### Refactored numeric integral implementation by Papagiannis et al.
```matlab
function d_Z = method_self_aeras_1_strwma(frequency, permittivity_layers, height)
    % Numerically evaluate self-impedance integrals using hybrid quadrature
    % frequency: frequency of operation [Hz]
    % permittivity_layers: permittivity values of the layers
    % height: distance from conductor to ground [m]

    % --- Constants ---
    tolerance = 1e-10;    % Convergence tolerance
    mu_0 = 4 * pi * 1e-7; % Magnetic permeability
    omega = 2 * pi * frequency;

    % Load Gauss-Legendre and Laguerre points/weights
    [s_legendre, w_legendre] = load_quadrature('ah16.txt');
    [s_laguerre, w_laguerre] = load_quadrature('ah35.txt');

    % --- Initialization ---
    u(1) = 4e-4 / (2 * height);
    dZ_previous = 0;
    iteration_limit = 100; % Fail-safe iteration cap

    % Compute initial segments
    Gu_1 = compute_legendre_integral(u(1), height, s_legendre, w_legendre, permittivity_layers, omega);
    Gu_2 = compute_laguerre_integral(u(1), height, s_laguerre, w_laguerre, permittivity_layers);
    dZ(1) = Gu_1 + Gu_2;

    % Iterative computation
    for i = 2:iteration_limit
        % Update step size
        step_factor = 1 / (2^(i - 1));
        u(i) = 10 * u(i - 1);

        % Gauss-Legendre segments
        Gu_3 = compute_legendre_integral(u(1) * step_factor, height, s_legendre, w_legendre, permittivity_layers, omega);
        Gu_4 = compute_legendre_integral(u(i) - u(i - 1), height, s_legendre, w_legendre, permittivity_layers, omega);

        % Laguerre contribution
        Gu_6 = compute_laguerre_integral(u(i), height, s_laguerre, w_laguerre, permittivity_layers);

        % Total impedance
        dZ(i) = Gu_3 + Gu_4 + Gu_6;

        % Check convergence
        if abs(dZ(i) - dZ_previous) <= tolerance
            break;
        end

        % Update for next iteration
        dZ_previous = dZ(i);
    end

    % Final result
    d_Z = 1j * omega * mu_0 * dZ(i) / pi;
end

% --- Helper Functions ---
function [points, weights] = load_quadrature(filename)
    % Load quadrature points and weights from file
    data = load(filename);
    points = data(:, 1);
    weights = data(:, 2);
end

function G = compute_legendre_integral(interval, height, points, weights, permittivity, omega)
    % Compute integral using Gauss-Legendre quadrature
    G = 0;
    for w = 1:length(weights)
        u_new = ((interval * (points(w) - 1)) / 2) + interval;
        F_u = nakagawa(u_new, omega, permittivity);
        G = G + weights(w) * exp(-2 * height * u_new) * F_u;
    end
    G = G * interval / 2;
end

function G = compute_laguerre_integral(u_start, height, points, weights, permittivity)
    % Compute integral using Laguerre quadrature
    G = 0;
    for v = 1:length(weights)
        u_new = points(v) / (2 * height) + u_start;
        F_u = nakagawa(u_new, 2 * pi * permittivity, permittivity);
        G = G + weights(v) * F_u;
    end
    G = G * exp(-2 * height * u_start) / (2 * height);
end

function F = nakagawa(u, omega, permittivity)
    % Compute the Nakagawa function
    mu_0 = 4 * pi * 1e-7;
    epsilon_0 = 8.854187817e-12;

    % Free-space propagation constant
    k_0 = omega * sqrt(mu_0 * epsilon_0);

    % Return Nakagawa function value
    F = 1 / (u + sqrt(u^2 + 1j * omega * mu_0 * (1 / permittivity)));
end
```


### **1. Initial computation using Gauss-Legendre quadrature**

#### Paper’s description:
> "The Gauss–Legendre method is applied to calculate the integral between zero and the first root... as this method is capable of handling the initial steep descent of the integrand."

#### Code mapping:
This is handled by the **first Gauss-Legendre integral** computation in the `method_self_aeras_1_strwma` function:

```matlab
Gu_1 = compute_legendre_integral(u(1), height, s_legendre, w_legendre, permittivity_layers, omega);
```

##### **`Gu_1`: Gauss-Legendre quadrature**
This term evaluates the integral over the interval $[0, u_1]$ using **Gauss-Legendre quadrature**, which is suited for **finite intervals**. The general description of the computation is as follows:
- Compute the initial segment of the integral from 0 to the first integration step size `u(1)`, where the function typically has a **steep descent** (as the paper mentions).
- **`compute_legendre_integral`** implements Gauss-Legendre quadrature, which uses roots of Legendre polynomials and weights for high accuracy in this segment.
- Gauss-Legendre quadrature is well-suited for this, as it accurately handles functions with rapid changes in this region.
- It divides $[0, u_1]$ into sub-intervals using the **Legendre points and weights**, which optimally sample the steep descent.

The first interval $u_1$ is estimated as:
```matlab
u(1) = 4e-4 / (2 * height);
```
This scales the step size by the conductor height, ensuring it is small enough to capture the steep descent behavior.

---

### **2. Evaluation using Laguerre quadrature**

#### Paper’s description:
> "Then the shifted Gauss-Laguerre method is used for the evaluation of the rest of the integral."

##### **`Gu_2`: Gauss-Laguerre Quadrature**
This term evaluates the **semi-infinite part of the integral** starting from $u_1$. Specifically, it uses **Gauss-Laguerre quadrature**, which is designed for intervals of the form $[u_1, \infty)$.

#### Code mapping:
This is handled by the **Laguerre quadrature** computation:

```matlab
Gu_2 = compute_laguerre_integral(u(1), height, s_laguerre, w_laguerre, permittivity_layers);
```

The general description of the computation is as follows:
- Compute the contribution of the integral for the semi-infinite part of the domain using Gauss-Laguerre quadrature, i.e. the part of the integral extending to **infinity**, starting at $u_1$.
- Laguerre quadrature naturally handles the **exponential decay** of terms like $e^{-2 h \cdot u}$, which dominate in the semi-infinite domain.
- Gauss-Laguerre is well-suited for semi-infinite intervals. The weights $w_\text{laguerre}$ and nodes $s_\text{laguerre}$ are designed to approximate integrals weighted by exponential decay terms $e^{-x}$.

The `compute_laguerre_integral` function implements:
```matlab
u_new = points(v) / (2 * height) + u_start;
```
- $u_{\text{start}} = u(1)$ ensures a shift to match the end of the Gauss-Legendre interval.
- $F_u$ is computed via the **Nakagawa function**, which is the main part of the integrand.

The nodes and weights of Laguerre quadrature are tailored for integrals involving **exponentials** like $e^{-x}$. The interval $[u_1, \infty)$ is effectively compressed into a manageable summation by Laguerre weights and nodes.

The paper explains why **two methods** are used:
1. **Gauss-Legendre for $[0, u_1]$**:
- This interval contains a **steep descent** of the integrand, which Gauss-Legendre handles efficiently for short, finite intervals.
- The sharp initial decay means that uniform quadrature (e.g., trapezoidal) would need many points to achieve similar accuracy.

2. **Gauss-Laguerre for $[u_1, \infty)$**:
- This region is dominated by the **exponential decay** $e^{-2 h \cdot u}$, making Laguerre quadrature the natural choice for efficient evaluation.
- Using Gauss-Legendre here would require transforming the infinite domain or truncating it, which is less efficient.

It is noted that **`Gu_1`** and **`Gu_2`** do not integrate over the same interval. They are complementary contributions to the overall integral:
- `Gu_1`: Covers $[0, u_1]$.
- `Gu_2`: Covers $[u_1, \infty)$.

#### **How Gauss-Laguerre handles semi-infinite domains**

Gauss-Laguerre quadrature is specifically designed for **semi-infinite integrals** of the form:

$\int_{0}^{\infty} e^{-x} f(x) \, dx$

To compute an integral over $[u_1, \infty)$, the **Laguerre transformation shifts the semi-infinite domain** $[0, \infty)$ to start from $u_1$. This shift is performed in the code using the term:

```matlab
u_new = points(v) / (2 * height) + u_start;
```

**Step 1: Original semi-infinite interval**
  The base **Laguerre quadrature nodes (`points`) and weights (`weights`)** are computed for the standard interval $[0, \infty)$, using the formula:

$\int_{0}^{\infty} e^{-u} f(u) \, du \approx \sum_{v=1}^n w_v f(x_v)$

Where:
- $x_v$: Gauss-Laguerre nodes.
- $w_v$: Gauss-Laguerre weights.

**Step 2: Shifting the interval to start at $u_1$**
  To extend the integration domain to $[u_1, \infty)$, the variable $u$ is **shifted** by adding $u_1$:

  $u_{\text{new}} = u_1 + \frac{x_v}{2 \cdot h}$

  From the code:
```matlab
u_new = points(v) / (2 * height) + u_start;
```
  
  - `u_start`: Acts as the starting point of the shifted interval $u_1$.
  - `points(v) / (2 * height)`: Scales the Laguerre nodes appropriately for the integration range.

  This transformation moves the integration domain from $[0, \infty)$ to $[u_1, \infty)$.

**Step 3: How the exponential decay is handled**

  Laguerre quadrature implicitly accounts for the exponential decay $e^{-u}$. When $u_{\text{new}} = u_1 + x_v / (2 \cdot h)$, the decay factor $e^{-2 h \cdot u_{\text{new}}}$ is naturally incorporated into the computation. The scaling factor $1 / (2 \cdot h)$ ensures proper handling of the semi-infinite range.
  
  From the code:
```matlab
G = G + weights(v) * F_u; % Sum the weighted function values
```

In the computation of `Gu_2`, the **Laguerre transformation shifts the effective interval** to $[u_1, \infty)$ without explicitly stating so. Instead, this happens "automatically" due to:

1. The use of shifted $u_{\text{new}} = u_1 + x_v / (2 \cdot h)$.
2. The exponential decay factor $e^{-u}$, inherent in Laguerre quadrature.

Thus, even though `Gu_2` is computed using Gauss-Laguerre points and weights, it effectively integrates from $u_1$ to infinity. This clever use of Gauss-Laguerre quadrature avoids having to directly handle the infinite upper bound $\infty$. The semi-infinite domain $[u_1, \infty)$ is "transformed" into a summation over Laguerre nodes and weights, with the scaling and shifting done via:

```matlab
u_new = points(v) / (2 * height) + u_start;
```

In conclusion, **`Gu_2` does indeed compute the contribution from $[u_1, \infty)$**, even though it’s not immediately obvious from the code.

---

### **3. Iterative extension of integration**

#### Paper’s description:
> "The procedure is repeated iteratively. In each iteration, the initial interval is bisected and the use of the Gauss–Legendre method is extended by intervals to the right of $u_0$."

#### Code mapping:
This is implemented in the **iterative loop** of `method_self_aeras_1_strwma`:

```matlab
for i = 2:iteration_limit
    % Update step size
    step_factor = 1 / (2^(i - 1));
    u(i) = 10 * u(i - 1);

    % Gauss-Legendre segments
    Gu_3 = compute_legendre_integral(u(1) * step_factor, height, s_legendre, w_legendre, permittivity_layers, omega);
    Gu_4 = compute_legendre_integral(u(i) - u(i - 1), height, s_legendre, w_legendre, permittivity_layers, omega);

    % Laguerre contribution
    Gu_6 = compute_laguerre_integral(u(i), height, s_laguerre, w_laguerre, permittivity_layers);

    % Total impedance
   

```matlab
    dZ(i) = Gu_3 + Gu_4 + Gu_6;

    % Check convergence
    if abs(dZ(i) - dZ_previous) <= tolerance
        break;
    end

    % Update for next iteration
    dZ_previous = dZ(i);
end
```

- **Purpose**:
  - Extends the integration by progressively adding contributions from new intervals, refining the solution iteratively.
  - The **step factor** halves the segment size at each iteration, improving resolution.
  - `u(i)` defines the new integration bound, increasing logarithmically to cover the semi-infinite domain.
- **Key steps**:
  1. **`Gu_3`**:
     - Computes a small segment at the leftmost interval $[0, u_1/2]$ using Gauss-Legendre.
  2. **`Gu_4`**:
     - Computes the interval between the two successive bounds $[u(i-1), u(i)]$ using Gauss-Legendre.
  3. **`Gu_6`**:
     - Adds the Laguerre contribution for the new interval $[u(i), \infty)$.
- **Convergence check**:
  - The code compares successive impedance estimates \( dZ(i) \) to \( dZ(i-1) \) using:
    ```matlab
    if abs(dZ(i) - dZ_previous) <= tolerance
    ```
    - The loop stops when the change falls below the predefined `tolerance`.

---

### **4. Final computation of impedance**

#### Paper’s description:
> "Convergence is achieved when the absolute difference between two successive values of the calculated integral is less than the predefined tolerance."

#### Code mapping:
This is reflected in the final line of the main function:

```matlab
d_Z = 1j * omega * mu_0 * dZ(i) / pi;
```

- **Purpose**:
  - Once convergence is achieved, the final self-impedance $\Delta Z$ is computed as:
    $\Delta Z = \frac{j \cdot \omega \cdot \mu_0 \cdot dZ(i)}{\pi}$
  - Here, `dZ(i)` represents the converged value of the integral.

---

### **Summary: how each part matches the paper**

| **Step in the paper**                                           | **Code implementation**                                                                                                                                   |
|------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Initial steep descent handled by Gauss-Legendre**             | `Gu_1 = compute_legendre_integral(u(1), ...)`                                                                                                             |
| **Shifted Gauss-Laguerre for semi-infinite domain**             | `Gu_2 = compute_laguerre_integral(u(1), ...)`                                                                                                             |
| **Iterative refinement via step bisection**                     | Iterative `for` loop, halving the step size: `step_factor = 1 / (2^(i-1))`, adding contributions with `Gu_3`, `Gu_4`, and `Gu_6`.                          |
| **Convergence check using tolerance**                           | `if abs(dZ(i) - dZ_previous) <= tolerance`, breaking the loop when the result stabilizes.                                                                 |
| **Final impedance computed after convergence**                  | `d_Z = 1j * omega * mu_0 * dZ(i) / pi`, assembling the result using the final converged value of the total integral.                                       |
