### **1. Initial computation using Gauss-Legendre quadrature**

#### Paper’s description:
> "The Gauss–Legendre method is applied to calculate the integral between zero and the first root... as this method is capable of handling the initial steep descent of the integrand."

#### Code mapping:
This is handled by the **first Gauss-Legendre integral** computation in the `method_self_aeras_1_strwma` function:

### **`Gu_1`: Gauss-Legendre quadrature**
This term evaluates the integral over the interval $[0, u_1]$ using **Gauss-Legendre quadrature**, which is suited for **finite intervals**.

```matlab
Gu_1 = compute_legendre_integral(u(1), height, s_legendre, w_legendre, permittivity_layers, omega);
```

- **Purpose**:
  - Compute the initial segment of the integral from 0 to the first integration step size \( u(1) \), where the function typically has a **steep descent** (as the paper mentions).
  - **`compute_legendre_integral`** implements Gauss-Legendre quadrature, which uses roots of Legendre polynomials and weights for high accuracy in this segment.
  - Gauss-Legendre quadrature is well-suited for this, as it accurately handles functions with rapid changes in this region.
    
- **Key behavior**:
  - It divides \( [0, u_1] \) into sub-intervals using the **Legendre points and weights**, which optimally sample the steep descent.
    
- **Details**:
  - \( u(1) \) is defined as:
    ```matlab
    u(1) = 4e-4 / (2 * height);
    ```
    This scales the step size by the height, ensuring it is small enough to capture the steep descent behavior.

---

### **2. Evaluation using Laguerre quadrature**

#### Paper’s description:
> "Then the shifted Gauss-Laguerre method is used for the evaluation of the rest of the integral."

### **`Gu_2`: Gauss-Laguerre Quadrature**
This term evaluates the **semi-infinite part of the integral** starting from \( u_1 \). Specifically, it uses **Gauss-Laguerre quadrature**, which is designed for intervals of the form \( [u_1, \infty) \).

#### Code mapping:
This is handled by the **Laguerre quadrature** computation:

```matlab
Gu_2 = compute_laguerre_integral(u(1), height, s_laguerre, w_laguerre, permittivity_layers);
```

- **Purpose**:
  - Compute the contribution of the integral for the semi-infinite part of the domain using Gauss-Laguerre quadrature.
  - Compute the part of the integral extending to **infinity**, starting at \( u_1 \).
  - Laguerre quadrature naturally handles the **exponential decay** of terms like \( e^{-2 h \cdot u} \), which dominate in the semi-infinite domain.

- **How it works**:
  - Gauss-Laguerre is well-suited for semi-infinite intervals. The weights \( w_\text{laguerre} \) and nodes \( s_\text{laguerre} \) are designed to approximate integrals weighted by exponential decay terms \( e^{-x} \).
  - The `compute_laguerre_integral` function implements:
    ```matlab
    u_new = points(v) / (2 * height) + u_start;
    ```
    - \( u_{\text{start}} = u(1) \) ensures a shift to match the end of the Gauss-Legendre interval.
    - \( F_u \) is computed via the **Nakagawa function**, which is the main part of the integrand.

  - **Key behavior**:
  - The nodes and weights of Laguerre quadrature are tailored for integrals involving **exponentials** like \( e^{-x} \). The interval \( [u_1, \infty) \) is effectively compressed into a manageable summation by Laguerre weights and nodes.

The paper clearly explains why **two methods** are used:
1. **Gauss-Legendre for \( [0, u_1] \)**:
   - This interval contains a **steep descent** of the integrand, which Gauss-Legendre handles efficiently for short, finite intervals.
   - The sharp initial decay means that uniform quadrature (e.g., trapezoidal) would need many points to achieve similar accuracy.

2. **Gauss-Laguerre for \( [u_1, \infty) \)**:
   - This region is dominated by the **exponential decay** \( e^{-2 h \cdot u} \), making Laguerre quadrature the natural choice for efficient evaluation.
   - Using Gauss-Legendre here would require transforming the infinite domain or truncating it, which is less efficient.

It is noted that **`Gu_1`** and **`Gu_2`** do not integrate over the same interval. They are complementary contributions to the overall integral:
- `Gu_1`: Covers \( [0, u_1] \).
- `Gu_2`: Covers \( [u_1, \infty) \).
  
---

### **3. Iterative Extension of Integration**

#### Paper’s Description:
> "The procedure is repeated iteratively. In each iteration, the initial interval is bisected and the use of the Gauss–Legendre method is extended by intervals to the right of \( u_0 \)."

#### Code Mapping:
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
- **Key Steps**:
  1. **`Gu_3`**:
     - Computes a small segment at the leftmost interval \( [0, u_1/2] \) using Gauss-Legendre.
  2. **`Gu_4`**:
     - Computes the interval between the two successive bounds \( [u(i-1), u(i)] \) using Gauss-Legendre.
  3. **`Gu_6`**:
     - Adds the Laguerre contribution for the new interval \( [u(i), \infty) \).
- **Convergence Check**:
  - The code compares successive impedance estimates \( dZ(i) \) to \( dZ(i-1) \) using:
    ```matlab
    if abs(dZ(i) - dZ_previous) <= tolerance
    ```
    - The loop stops when the change falls below the predefined `tolerance`.

---

### **4. Final Computation of Impedance**

#### Paper’s Description:
> "Convergence is achieved when the absolute difference between two successive values of the calculated integral is less than the predefined tolerance."

#### Code Mapping:
This is reflected in the final line of the main function:

```matlab
d_Z = 1j * omega * mu_0 * dZ(i) / pi;
```

- **Purpose**:
  - Once convergence is achieved, the final self-impedance \( \Delta Z \) is computed as:
    \[
    \Delta Z = \frac{j \cdot \omega \cdot \mu_0 \cdot dZ(i)}{\pi}
    \]
  - Here, \( dZ(i) \) represents the converged value of the integral.

---

### **Summary: How Each Part Matches the Paper**

| **Step in the Paper**                                           | **Code Implementation**                                                                                                                                   |
|------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Initial steep descent handled by Gauss-Legendre**             | `Gu_1 = compute_legendre_integral(u(1), ...)`                                                                                                             |
| **Shifted Gauss-Laguerre for semi-infinite domain**             | `Gu_2 = compute_laguerre_integral(u(1), ...)`                                                                                                             |
| **Iterative refinement via step bisection**                     | Iterative `for` loop, halving the step size: `step_factor = 1 / (2^(i-1))`, adding contributions with `Gu_3`, `Gu_4`, and `Gu_6`.                          |
| **Convergence check using tolerance**                           | `if abs(dZ(i) - dZ_previous) <= tolerance`, breaking the loop when the result stabilizes.                                                                 |
| **Final impedance computed after convergence**                  | `d_Z = 1j * omega * mu_0 * dZ(i) / pi`, assembling the result using the final converged value of the total integral.                                       |
