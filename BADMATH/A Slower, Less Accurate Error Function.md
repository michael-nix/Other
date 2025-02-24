# A Slower, Less Accurate Error Function
Doing probability theory is too easy, let's make it worse.

**Michael Nix**, Toronto, Canada, 2025

---

The only thing I remember from engineering school is every single Fourier Transform pair.  As though etched into my bones, I can't actually do the math, only grasp around and guess at what ends up being something pretty close to the right answer.  I mean, skinny functions transform to wide ones, finite sharp functions become smooth and infinite, or like, squares go to sincs.  Convolution or multiplication of any one of these by another moves them around or reshapes them.  Ya know, the basics.  Which is how this came to me in a dream:

```matlab
gx = @(x, n, s, T) 2/T * (0.5 + sum(exp(-(2*pi/T*s * (1:n).').^2/2) .* cos((1:n).' .* (2*pi/T * x)), 1))
```

It turns out that one of my hobbies is dreaming up terrible one line functions in MATLAB that mimic real math.  Or rather, make worse versions of built in math functions.  Slower, less accurate.  This one is a normalized Gaussian function made up from a sum of discrete Fourier series components.  Since Gaussians are generally symmetric about zero, $n$ is the number of cosines used to construct it, $s$ is $\sigma$, the standard deviation, and $T$ is the period of the series (since Fourier series repeat).  

Other than the deepest recesses of my mind, where in the hell did this come from?

When you look up, "how to do a Fourier transform for a Gaussian," on the internet, the first thing you see is the Wolfram MathWorld answer which boils down to, "it's just a cosine transfer because sines are odd, and then you look it up in Abramowitz and Stegun."  This is not a satisfying answer.  There are some standard derivations out there, working through the actual math nicely, but they're boring and I found them hard to follow.  Which is why it took me so long to type this up.  But then I found a good one posted online by Konstantinos Derpanis at York University.  Here it is:

First, we're going to use this definition of the Fourier transform,

```math
F(\omega) = \int_{-\infty}^\infty f(t) \, \mathrm{e}^{-j\omega t} \mathrm{d} t
```

Then, starting with the definition of a normalized Gaussian,

```math
f(t) = \frac{1}{\sqrt{2\pi\sigma^2}} \mathrm{e}^{-x^2/2\sigma^2}
```

We take the derivative of both sides,

```math
\begin{aligned}
\frac{\mathrm{d}}{\mathrm{d}t} f(t) &= -\frac{x}{\sigma^2} \frac{1}{\sqrt{2\pi\sigma^2}} \mathrm{e}^{-x^2/2\sigma^2} \\
&= -\frac{t}{\sigma^2} f(t)
\end{aligned}
```

Which is always fun, since we always get the original Gaussian back (this is where Hermite polynomials come from).  We then take the Fourier transform of both sides,

```math
\begin{aligned}
j\omega F(\omega) &= -\frac{t}{\sigma^2} F(\omega) \\
&= -\frac{j}{\sigma^2} \frac{\mathrm{d}}{\mathrm{d}\omega} F(\omega)
\end{aligned}
```

And re-arrange so that all of our Fourier transformed parts are on the left-hand side,

```math
\frac{1}{F(\omega)} \frac{\mathrm{d}}{\mathrm{d}\omega} F(\omega) = -\omega\sigma^2
```

We then note, that via the chain rule, the derivative of the logarithm of a Fourier transform is the same thing as what we already have on the left-hand side, that is,

```math
\frac{\mathrm{d}}{\mathrm{d}\omega} \log F(\omega) = \frac{1}{F(\omega)} \frac{\mathrm{d}}{\mathrm{d}\omega} F(\omega)
```

So we make the substitution back, giving us,

```math
\frac{\mathrm{d}}{\mathrm{d}\omega} \log F(\omega) = -\omega\sigma^2
```

We can then integrate with respect to $\omega$ using the fundamental theorem of calculus,

```math
\int_0^\omega\frac{\mathrm{d}}{\mathrm{d}\omega'} \log F(\omega') \;\mathrm{d}\omega' = -\sigma^2 \int_0^\omega \omega' \;\mathrm{d}\omega'
```

Which leads us to the simple expression,

```math
\log F(\omega) - \log F(0) = -\frac{1}{2}\sigma^2\omega^2
```

But, $F(0)$ is just the Fourier transform of our function when $\omega$ is set to zero.  So using our definition of the Fourier transform, the exponential of the $\omega$ terms become one, meaning we're just integrating our Gaussian,

```math
F(0) = \int_{-\infty}^\infty f(t) \; \mathrm{d} t = 1
```

Which of course is already properly normalized, and the logarithm of one is just zero, so we're left with,

```math
\log F(\omega) = -\frac{1}{2}\sigma^2\omega^2
```

We can then raise $\mathrm{e}$ to the power of both sides, giving us,

```math
\mathrm{e} ^{\log F(\omega)} = \mathrm{e}^{-\frac{1}{2}\sigma^2\omega^2}
```

And since $\mathrm{e}$ raised to the power of the logarithm of some value is just that value, we're left with our final answer:

```math
\boxed{\therefore \; F(\omega) = \mathrm{e}^{-\frac{1}{2}\sigma^2\omega^2}}
```

Which is pretty clever!  The hardest integration we had to do was to remember that the integral of a special function was one, nice.  But we're only halfway there, this is the Fourier transform of a Gaussian, but how do we use this to construct Gaussians from cosines?  

The first step is to use the definition of the inverse transform, realize that sines are again odd functions, and then convert from angular frequency to regular frequency.  This gets us,

```math
\begin{aligned}
f(t) &= \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2\omega^2} \mathrm{e}^{j\omega t} \; \mathrm{d}\omega \\
&= \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2\omega^2} \cos(\omega t) \; \mathrm{d}\omega \\
&= \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi f)^2} \cos(2\pi f t) \; \mathrm{d}f
\end{aligned}
```

If we then convert this into a Fourier _series_, it will then become:

```math
f(t) \approx \frac{1}{T} \sum_{n\,=\,-\infty}^{\infty} \mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \cos\left(2\pi \frac{n}{T} t\right)
```

Which looks like this, when dealing with $\mathcal{N}(0, 0.25)$:

<p align="center">
    <img src="./figures/periodic gaussian.png"><br>
    <i>Figure 1: Periodic Gaussian by Fourier Series Approximation</i>
</p>

And I hope you notice that--much like every Fourier series transform--this is a periodic approximation.  That is, the Gaussian we have here is repeated for every interval, $T$.  

But this isn't what we set out to do.  We wanted to create a bad version of a perfectly well approximated function, the error function.  In order to get us closer to our true goal, we have to do two things:

1. Correct it so that it's no longer periodic,
2. Integrate it to get a cumulative distribution function,
3. Translate that to an error function.

```math
-\frac{T}{2} \le t \le \frac{T}{2}
```

```math
T = 2k\,\sigma
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
\Phi \left ( \frac{t-\mu}{\sigma} \right ) \approx \frac{2}{T} \left (\frac{t-\mu}{2\sigma} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} (2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T\sigma} (t-\mu)\right) \right ) + \frac{1}{2}
```

```math
\Phi \left ( \frac{t-\mu}{\sigma} \right ) = \frac{1}{2} \left [1 + \mathrm{erf} \left ( \frac{t-\mu}{\sigma \sqrt{2}} \right ) \right ]
```

```math
\begin{aligned}
\mathrm{erf} \left ( \frac{t-\mu}{\sigma \sqrt{2}} \right ) &= 2 \Phi \left ( \frac{t-\mu}{\sigma} \right ) - 1 \\

&\approx \frac{4}{T} \left (\frac{t-\mu}{2\sigma} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} (2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T\sigma} (t-\mu)\right) \right )
\end{aligned}
```