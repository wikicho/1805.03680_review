# 1805.03680_review

Bounce action is also can be calculated with auxilary function $V_t(\phi)$~\cite{Espinosa:2018hue}.
In terms of $V$ and $V_t$ this action can be written as
$$
\begin{equation}
    S_{E,d}=\frac{(d-1)^{(d-1)}(2\pi)^{d/2}}{\Gamma(1+d/2)}\int_{\phi_+}^{\phi_0}\frac{(V-V_t)^{d/2}}{\left|V_t'\right|^{d-1}}d\phi\ .
\end{equation}
$$


In $d$-dimension (where $d = 3,4$), linear one is 
$$
\begin{equation}
V_{t1}(\phi)=V_0\frac{\phi}{\phi_0}\ ,
\end{equation}
$$
where $V_0=V(\phi_0)$, and it satisfies the boundary conditions (\eqref{eq:chap5:vt_boundary_condition}). 
To also satisfy the boundary condition on $V_t'(\phi_0)$, 
$$
\begin{equation}
V_t'(\phi_+) = 0\ ,\quad V_t'(\phi_0) = \frac{(d-1)}{d} V'(\phi_0)\ .
\label{eq:chap5:vt_boundary_condition}
\end{equation}
$$
one can use the quadratic approximation
$$
\begin{equation}
V_{t2}(\phi)=V_{t1}(\phi)+\frac{\phi}{d\phi_0^2}\left[(d-1)\phi_0V'_0-dV_0\right](\phi-\phi_0)\ ,
\end{equation}
$$
where $V'_0=V'(\phi_0)$. One can go further, matching also $V_t'(\phi_+)$ using the cubic approximation
$$
\begin{equation}
V_{t3}(\phi)=V_{t2}(\phi)+\frac{\phi}{d\phi_0^3}\left[(d-1)\phi_0 V'_0-2d V_0\right](\phi-\phi_0)^2\ .
\end{equation}
$$
where we used the boundary condition and $\phi_{\text{FV}} = 0$. Adding the information about the barrier that separating the vacua, $\phi = \phi_T$ where local maximum of effective potential in range of $(\phi_{FV}, \phi_0)$. It can be achieved with the quartic approximation
$$
\begin{equation}
    V_{t4}(\phi) = V_{t3}(\phi) + q_d \phi^{2} (\phi-\phi_0)^{2}.
\end{equation}
$$
Start with \eqref{eq:chap5:vt_boundary_condition}
$$
\begin{equation}
d(V_{t}^{\prime})^2  = 2(d-1)\bigl[V-V_t\bigr]V_t^{\prime\prime} + (d-1) V^{\prime} V_t^{\prime}
\label{eq:EoMd}
\end{equation}
$$
Setting 
\begin{align*}
 Q(\phi) &\equiv \phi^{2} (\phi-\phi_0)^{2}.
\end{align*}
$$
At $\phi=\phi_T$, equation~\eqref{eq:EoMd} reduces to
$$
\begin{equation}
d (V_{t3}'(\phi_T) + q_d Q_{T}^{\prime})^2 + 
 2 (d - 1) (V(\phi_T) - V_{t3}(\phi_T) - q_d Q_T) (V_{t3}''(\phi_T)+ q_d Q_{T}^{\prime\prime}) = 0 .
\end{equation}
$$
Collecting powers of $q_d$ one finds a quadratic
$$
\begin{equation}
a_d q_d^{2} + 2 b_d q_d + c_d = 0,
\end{equation}
$$
where
$$
\begin{equation}
    \begin{split}
        a_d &= d {Q^{\prime}_T}^{2} - 2(d-1) Q_T Q_T^{\prime\prime}, \\
        b_d &= -2 C (d-1) Q_T + 2 d B  Q_T^{\prime} - 2 (d - 1) Q_T^{\prime\prime} (A - V(\phi_T)), \\
        c_d &= B^2  - 2 C (d-1) (A - V(\phi_T))
    \end{split}
\end{equation} 
$$
Choosing the branch that continuously connects to $q_d\to0$ when $V_{t3}$ is exact,
\begin{equation}
q_d = \frac{1}{a_d} \left(-b_d - \sqrt{b_d^2 - c_d a_d} \right).
\label{eq:a4d}
\end{equation}
$$
