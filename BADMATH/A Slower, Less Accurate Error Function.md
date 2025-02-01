# A Slower, Less Accurate Error Function
Doing probability theory is too easy, let's make it worse.

**Michael Nix**, Toronto, Canada, 2025

---



```matlab
gx = @(x, n, s, T) 2/T * (0.5 + sum(exp(-(2*pi/T*s * (1:n).').^2/2) .* cos((1:n).' .* (2*pi/T * x)), 1))
```

```math
F(\omega) = \int_{-\infty}^\infty f(t) \, \mathrm{e}^{-j\omega t} \mathrm{d} t
```

```math
f(t) = \frac{1}{\sqrt{2\pi\sigma^2}} \mathrm{e}^{-x^2/2\sigma^2}
```

```math
\begin{aligned}
\frac{\mathrm{d}}{\mathrm{d}t} f(t) &= -\frac{x}{\sigma^2} \frac{1}{\sqrt{2\pi\sigma^2}} \mathrm{e}^{-x^2/2\sigma^2} \\
&= -\frac{t}{\sigma^2} f(t)
\end{aligned}
```

```math
\begin{aligned}
j\omega F(\omega) &= -\frac{t}{\sigma^2} F(\omega) \\
&= -\frac{j}{\sigma^2} \frac{\mathrm{d}}{\mathrm{d}\omega} F(\omega)
\end{aligned}
```

```math
\frac{1}{F(\omega)} \frac{\mathrm{d}}{\mathrm{d}\omega} F(\omega) = -\omega\sigma^2
```

```math
\frac{\mathrm{d}}{\mathrm{d}\omega} \log F(\omega) = \frac{1}{F(\omega)} \frac{\mathrm{d}}{\mathrm{d}\omega} F(\omega)
```

```math
\frac{\mathrm{d}}{\mathrm{d}\omega} \log F(\omega) = -\omega\sigma^2
```

```math
\int_0^\omega\frac{\mathrm{d}}{\mathrm{d}\omega'} \log F(\omega') \;\mathrm{d}\omega' = -\sigma^2 \int_0^\omega \omega' \;\mathrm{d}\omega'
```

```math
\log F(\omega) - \log F(0) = -\frac{1}{2}\sigma^2\omega^2
```

```math
F(0) = \int_{-\infty}^\infty f(t) \; \mathrm{d} t = 1
```

```math
\log F(\omega) = -\frac{1}{2}\sigma^2\omega^2
```

```math
\mathrm{e} ^{\log F(\omega)} = \mathrm{e}^{-\frac{1}{2}\sigma^2\omega^2}
```

```math
\boxed{\therefore \; F(\omega) = \mathrm{e}^{-\frac{1}{2}\sigma^2\omega^2}}
```

```math
\begin{aligned}
f(t) &= \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2\omega^2} \mathrm{e}^{j\omega t} \; \mathrm{d}\omega \\
&= \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2\omega^2} \cos(\omega t) \; \mathrm{d}\omega \\
&= \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi f)^2} \cos(2\pi f t) \; \mathrm{d}f
\end{aligned}
```

```math
f(t) \approx \frac{1}{T} \sum_{n\,=\,-\infty}^{\infty} \mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \cos\left(2\pi \frac{n}{T} t\right)
```

```math
-\frac{T}{2} \le t \le \frac{T}{2}
```

```math
T = k\,\sigma
```

```math
f(t) \approx \frac{2}{T} \left (\frac{1}{2} + \sum_{n\,=\,1}^{n_\mathrm{max}} \mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \cos\left(2\pi \frac{n}{T} t\right) \right )
```

```math
\Phi(t) = \int f(t) \; \mathrm{d}t \approx \frac{2}{T} \left (\frac{t}{2} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T} t\right) \right ) + C
```

```math
\Phi(t) \approx \frac{2}{T} \left (\frac{t}{2} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T} t\right) \right ) + \frac{1}{2}
```

```math
\begin{aligned}
f(t - \mu) &= \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2\omega^2} \mathrm{e}^{-j\omega \mu} \mathrm{e}^{j\omega t} \; \mathrm{d}\omega \\

&= \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2\omega^2} \mathrm{e}^{j\omega (t-\mu)} \; \mathrm{d}\omega \\

&= \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi f)^2} \cos(2\pi f (t-\mu)) \; \mathrm{d}f
\end{aligned}
```

```math
\Phi \left ( \frac{t-\mu}{\sigma} \right ) \approx \frac{2}{T} \left (\frac{t-\mu}{2} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T} (t-\mu)\right) \right ) + \frac{1}{2}
```

```math
\Phi \left ( \frac{t-\mu}{\sigma} \right ) = \frac{1}{2} \left [1 + \mathrm{erf} \left ( \frac{t-\mu}{\sigma \sqrt{2}} \right ) \right ]
```

```math
\begin{aligned}
\mathrm{erf} \left ( \frac{t-\mu}{\sigma \sqrt{2}} \right ) &= 2 \Phi \left ( \frac{t-\mu}{\sigma} \right ) - 1 \\

&\approx \frac{4}{T} \left (\frac{t-\mu}{2} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T} (t-\mu)\right) \right )
\end{aligned}
```