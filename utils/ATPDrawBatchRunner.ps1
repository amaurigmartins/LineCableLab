# =========================================================================
# ATPDraw Batch Runner
#
# DESCRIPTION:
#   This script automates running ATPDraw simulations for all specified
#   files within a directory and its subdirectories. It opens each file,
#   simulates an 'F2' keypress to start the simulation, monitors the output
#   .pl4 file until it's complete, and then closes ATPDraw before
#   processing the next file.
#
# HOW TO USE:
#   1. Configure the paths and parameters in the "MAIN SETTINGS" section.
#   2. Run the script from a PowerShell terminal.
#   3. Press 'P' during execution to pause/resume the script.
#
# =========================================================================

# --- MAIN SETTINGS ---

# 1. Paths for Executables and Files
$exePath = "C:\ATPDraw\Atpdraw.exe"
$workingDirectory = "C:\arq_ATP\EMI_protection"
# File patterns to search for recursively (can be a single pattern or an array)
$filePatterns = @("*.xml", "*.acp")

# 2. ATP Monitoring & Timeout Parameters
#    Directory where ATPDraw creates its temporary output files.
$atpWorkDir = "C:\ATP\work"
# The target file size in bytes that indicates a simulation is complete.
$targetPl4SizeInBytes = 21095e3
# How many seconds to wait AFTER the target size is reached before closing ATPDraw.
$postSizeWaitTime = 2
# How many seconds to wait for the .pl4 file to appear after pressing F2.
$waitTimeForFileCreation = 5
# Timeout in seconds if the .pl4 file never reaches the target size.
$maxWaitTimeForSize = 60
# How many times to retry pressing F2 if the .pl4 file doesn't appear.
$maxF2Retries = 3

# 3. Logging
#    A log file will be created here listing any simulations that failed.
$failedFilesLogPath = Join-Path $workingDirectory "Files_failed.txt"

# --- END OF SETTINGS ---


# --- SCRIPT INITIALIZATION ---

# Delete old log file if it exists to ensure a clean run
if (Test-Path $failedFilesLogPath) { Remove-Item $failedFilesLogPath }

# Define a C# class within PowerShell for advanced UI control (like bringing a window to the front)
if (-not ([System.Management.Automation.PSTypeName]'WinAPI').Type) {
    Add-Type -TypeDefinition @"
    using System;
    using System.Runtime.InteropServices;
    public class WinAPI {
        [DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr hWnd);
        [DllImport("user32.dll")] public static extern short GetAsyncKeyState(int vKey);
        public const int VK_P = 0x50; // Virtual-Key code for 'P'
    }
"@
}
# Load the necessary assembly for sending keystrokes
Add-Type -AssemblyName System.Windows.Forms

# --- PAUSE/RESUME FUNCTIONALITY ---
# This function allows you to pause the script by pressing the 'P' key.
$script:isPaused = $false
$script:pKeyWasDown = $false
function Check-Pause {
    # Check if the 'P' key is currently pressed
    $isKeyDown = ([WinAPI]::GetAsyncKeyState([WinAPI]::VK_P) -lt 0)
    # If it was just pressed (and wasn't down before), toggle the pause state
    if ($isKeyDown -and -not $script:pKeyWasDown) {
        $script:isPaused = -not $script:isPaused
        if ($script:isPaused) { Write-Host "`n--- SCRIPT PAUSED ---`nPress 'P' again to resume." -ForegroundColor Magenta }
        else { Write-Host "`n--- SCRIPT RESUMED ---`n" -ForegroundColor Magenta }
    }
    $script:pKeyWasDown = $isKeyDown

    # Loop here as long as the script is paused
    while ($script:isPaused) {
        $isKeyDownNow = ([WinAPI]::GetAsyncKeyState([WinAPI]::VK_P) -lt 0)
        if ($isKeyDownNow -and -not $script:pKeyWasDown) {
            $script:isPaused = $false
            Write-Host "`n--- SCRIPT RESUMED ---`n" -ForegroundColor Magenta
        }
        $script:pKeyWasDown = $isKeyDownNow
        Start-Sleep -Milliseconds 100
    }
}


# --- MAIN EXECUTION LOOP ---

# 1. Find all files matching the patterns recursively
Write-Host "Searching for files in '$workingDirectory'..." -ForegroundColor Cyan
$filesToProcess = Get-ChildItem -Path $workingDirectory -Include $filePatterns -Recurse
$totalFiles = $filesToProcess.Count

if ($totalFiles -eq 0) {
    Write-Host "No files found matching the specified patterns. Exiting." -ForegroundColor Yellow
    exit
}

Write-Host "Found $totalFiles files to process." -ForegroundColor Green
Start-Sleep -Seconds 2

# Initialize tracking variables
$stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
$failedSimulations = [System.Collections.ArrayList]@()
$currentIteration = 0

# 2. Loop through each file found
foreach ($file in $filesToProcess) {
    $currentIteration++
    Check-Pause
    Write-Host "--- Starting file $currentIteration of $totalFiles: $($file.FullName) ---" -ForegroundColor Yellow

    # Define the expected output file path
    $baseName = [System.IO.Path]::GetFileNameWithoutExtension($file.FullName)
    $expectedPl4File = Join-Path $atpWorkDir "$baseName.pl4"

    # Clean up any old output file from a previous run
    if (Test-Path $expectedPl4File) {
        Write-Host "[$currentIteration] Deleting old .pl4 file..." -ForegroundColor DarkGray
        Remove-Item $expectedPl4File -Force
    }

    # 3. Open the file in ATPDraw
    Write-Host "[$currentIteration] Opening file in ATPDraw..."
    Start-Process -FilePath $exePath -ArgumentList "`"$($file.FullName)`""
    Start-Sleep -Seconds 8 # Wait for the application to load

    # 4. Find the ATPDraw process and send the 'F2' keypress
    $proc = Get-Process | Where-Object { $_.MainWindowTitle -like "*ATPDraw*" } | Select-Object -First 1
    if (-not $proc) {
        Write-Host "[$currentIteration] ERROR: Could not find the ATPDraw process. Skipping." -ForegroundColor Red
        [void]$failedSimulations.Add($file.FullName)
        continue
    }

    [void][WinAPI]::SetForegroundWindow($proc.MainWindowHandle)
    Start-Sleep -Milliseconds 500

    # 5. Attempt to run simulation and monitor for the .pl4 file
    $pl4Created = $false
    for ($attempt = 1; $attempt -le $maxF2Retries; $attempt++) {
        Check-Pause
        Write-Host "[$currentIteration] Sending F2 to start simulation (Attempt $attempt of $maxF2Retries)..."
        [System.Windows.Forms.SendKeys]::SendWait("{F2}")
        Write-Host "[$currentIteration] Waiting for .pl4 file creation..."
        Start-Sleep -Seconds $waitTimeForFileCreation

        if (Test-Path $expectedPl4File) {
            Write-Host "[$currentIteration] .pl4 file detected successfully!" -ForegroundColor Green
            $pl4Created = $true
            break # Exit the retry loop
        } else {
            Write-Host "[$currentIteration] WARNING: .pl4 file not detected after attempt $attempt." -ForegroundColor Yellow
        }
    }

    if (-not $pl4Created) {
        Write-Host "[$currentIteration] ERROR: .pl4 file was not created after $maxF2Retries attempts. Skipping." -ForegroundColor Red
        [void]$failedSimulations.Add($file.FullName)
        Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
        continue
    }

    # 6. Wait for the simulation to complete by monitoring file size
    Write-Host "[$currentIteration] Monitoring .pl4 file until it reaches target size (Timeout: $maxWaitTimeForSize s)..."
    $sizeWatch = [System.Diagnostics.Stopwatch]::StartNew()
    $sizeReached = $false
    while ($sizeWatch.Elapsed.TotalSeconds -lt $maxWaitTimeForSize) {
        Check-Pause
        $currentSize = (Get-Item $expectedPl4File -ErrorAction SilentlyContinue).Length
        if ($currentSize -ge $targetPl4SizeInBytes) {
            Write-Host "[$currentIteration] Target file size reached." -ForegroundColor Green
            $sizeReached = $true
            break
        }
        Start-Sleep -Milliseconds 500 # Check twice per second
    }
    $sizeWatch.Stop()

    if (-not $sizeReached) {
        Write-Host "[$currentIteration] TIMEOUT: .pl4 file did not reach target size. Skipping." -ForegroundColor Red
        [void]$failedSimulations.Add($file.FullName)
        Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
        continue
    }

    # Wait a moment after completion for file handles to be released
    Write-Host "[$currentIteration] Waiting for post-simulation delay..."
    Start-Sleep -Seconds $postSizeWaitTime

    # 7. Close ATPDraw
    Write-Host "[$currentIteration] Closing ATPDraw..."
    Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
    Start-Sleep -Seconds 2 # Wait for the process to close fully
}


# --- FINAL SUMMARY ---
$stopwatch.Stop()
$elapsedTime = $stopwatch.Elapsed
$formattedTime = "{0:hh\:mm\:ss}" -f $elapsedTime

Write-Host ""
Write-Host "--- AUTOMATION PROCESS COMPLETED ---" -ForegroundColor Green
Write-Host "Total execution time: $formattedTime"

# Save the list of failed simulations to a text file, if any
if ($failedSimulations.Count -gt 0) {
    $logHeader = "The following $($failedSimulations.Count) simulations failed:"
    $logHeader | Set-Content -Path $failedFilesLogPath
    $failedSimulations | ForEach-Object { Add-Content -Path $failedFilesLogPath -Value "- $_" }

    Write-Host ""
    Write-Host "--- FAILED SIMULATIONS REPORT ---" -ForegroundColor Red
    Write-Host "A list of failed simulations has been saved to: $failedFilesLogPath"
} else {
    Write-Host "All simulations completed successfully." -ForegroundColor Green
}
