# A Slower, Less Accurate Error Function
Doing probability theory is too easy, let's make it worse.

**Michael Nix**, Toronto, Canada, 2025

---

The only thing I remember from engineering school is every single Fourier Transform pair.  As though etched into my bones, I can't actually do the math, only grasp around and guess at what ends up being something pretty close to the right answer.  I mean, skinny functions transform to wide ones, finite sharp functions become smooth and infinite, or like, squares go to sincs.  Convolution or multiplication of any one of these by another moves them around or reshapes them.  Ya know, the basics.  Which is how this came to me in a dream:

```matlab
gx = @(x, n, s, T) 2/T * (0.5 + sum(exp(-(2*pi/T*s * (1:n).').^2/2) .* cos((1:n).' .* (2*pi/T * x))))
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
    <img src="./figures/periodic gaussian.png" width="75%"><br>
    <i>Figure 1: Periodic Gaussian by Fourier Series Approximation</i>
</p>

And I hope you notice that--much like every Fourier series transform--this is a periodic approximation.  That is, the Gaussian we have here is repeated for every interval, $T$.

But this isn't what we set out to do.  We wanted to create a bad version of a perfectly well approximated function, the error function.  In order to get us closer to our true goal, we have to do two things:

1. Correct it so that it's no longer periodic,
2. Integrate it to get a cumulative distribution function,
3. Translate that to an error function.

The first thing we have to do is get rid of that periodicity.  Which is very easy when using a computer: just ignore it!  That is, if we limit the domain of our function to be within $\pm T / 2$, where $T$ is the period of our approximation, i.e.:

```math
-\frac{T}{2} \le t \le \frac{T}{2}
```

Then all we have to do is choose a suitable period.  Say we want some kind of use for this out to five standard deviations, then all we have to do in order to make sure we get both sides of the function working, is to double it:

```math
T = 2k\,\sigma
```

And in this way we make sure that we have $k$ standard deviations on either side of our function.  Before we go any further, we can first make a few simplifcations: 1) realize that the cosine of zero is just one; and, 2) realize that cosines are even functions, so we double their value, but only count from $n=1$ to infinity, and correct for double-counting of the $n=0$ case:

```math
f(t) \approx \frac{2}{T} \left (\frac{1}{2} + \sum_{n\,=\,1}^{n_\mathrm{max}} \mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \cos\left(2\pi \frac{n}{T} t\right) \right )
```

This is where that first MATLAB one-liner comes from.

Now we're only a few steps away from our final goal: some kind of error function.  The nice thing about our current approximation is that it's incredibly linear.  That is, adding up a bunch of cosines is almost the definition of linear--not at all, but still.  However, this means that we can do one thing with this approximation that we can't really do with normal Gaussians: just integrate the fucker.  Yup, that's it, we just turn the cosines into sines and shuffle the arguments around together, and this gives us:

```math
\Phi(t) = \int f(t) \; \mathrm{d}t \approx \frac{2}{T} \left (\frac{t}{2} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T} t\right) \right ) + C
```

But uh oh, that constant of integration has come to taunt us.  But thankfully cumulative distribution functions have really easy boundary conditions: at left infinity they're zero, and at right infinity they're one.  So we can just take a look at a graph of our function, and see where the ends end up.

<p align="center">
    <img src="./figures/cdf offset determination.png" width="75%"><br>
    <i>Figure 2: Figuring Out Constants of Integration</i>
</p>

Ok, because we summed up a bunch of sines, it's no surprise that at $t=0$ our approximation is going to be zero, but our ends don't go to any value.  Because we integrated a periodic function, the cumulative distribution function doesn't level off at either end, but at left infinity it tends to bottom infinity, and at right infinity it tends to top infinity.  That's not great.  But, in this graph we obviously didn't do what we said we were going to do earlier: we didn't truncate the domain of the function so that it's well defined in our region of interest.  So, if we squint and look at, say, $\pm$ 5, we can tell that our constant of integration is obviously $\frac{1}{2}$.  So that leaves us with a final form for our Gaussian CDF approximation, which we will call `badcdf` in this folder:

```math
\Phi(t) \approx \frac{2}{T} \left (\frac{t}{2} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T} t\right) \right ) + \frac{1}{2}
```

Or as a MATLAB one-liner:

```matlab
@(x, n, s, T) 2/T * (x/2 + sum(exp(-s^2 * ((1:n).').^2/2 * (2*pi/T)^2) .* sin((1:n).' .* x * (2*pi/T)) ./ ((1:n).' * (2*pi/T)))) + 1/2;
```

Where `x` is our input parameter, `n` is the number of sines to use, `s` is the standard deviation, and `T` is the period.

One of the interesting things about the error function is that it's usually given in terms of some kind of location and scale parameter (like a mean and standard deviation), not as separate inputs to the function, but as direct modifications to your input.  This means, that if we set our previous standard deviation parameter, $\sigma$ to one, we can reframe everything similar, meaning our cumulative distribution function then becomes:

```math
\Phi \left ( \frac{t-\mu}{\sigma} \right ) \approx \frac{2}{T} \left (\frac{t-\mu}{2\sigma} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} (2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T\sigma} (t-\mu)\right) \right ) + \frac{1}{2}
```

Now, this wouldn't work with our previous Gaussian apprxoimation, as if we set $\sigma$ to one, then the normalization gets all messed up.  But when we integrated it, the shuffling around of parameters actually took care of that for us.  This now means that we can directly compare our CDF approximation with the definition of a Gaussian CDF using the error function:

```math
\Phi \left ( \frac{t-\mu}{\sigma} \right ) = \frac{1}{2} \left [1 + \mathrm{erf} \left ( \frac{t-\mu}{\sigma \sqrt{2}} \right ) \right ]
```

And then plot them to get some indication that they're approximately the same:

<p align="center">
    <img src="./figures/cdf comparison.png" width="75%"><br>
    <i>Figure 3: Comparison of Gaussian CDF and Approximation</i>
</p>

As you can see, just directly scaling and relocating the input parameter for both functions gets what you'd expect for a Gaussian CDF with a mean of 1, and standard deviation of half.

Convinced that we're on the right track, we can now finally do what we've set out to do: define our slower, less accurate error function, included as `baderf` somewhere in the folder that this document is placed:

<!-- ```math
\begin{aligned}
\mathrm{erf} \left ( \frac{t-\mu}{\sigma \sqrt{2}} \right ) &= 2 \Phi \left ( \frac{t-\mu}{\sigma} \right ) - 1 \\
&\approx \frac{4}{T} \left (\frac{t-\mu}{2\sigma} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} (2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T\sigma} (t-\mu)\right) \right )
\end{aligned}
``` -->

```math
\begin{aligned}
\mathrm{erf}  ( t ) &= 2 \Phi \left ( t\sqrt{2} \right ) - 1 \\
&\approx \frac{4}{T} \left (\frac{t\sqrt{2}}{2} + \sum_{n\,=\,1}^{n_\mathrm{max}} \frac{T}{2\pi n}\mathrm{e}^{-\frac{1}{2} (2\pi \frac{n}{T})^2} \sin\left(2\pi \frac{n}{T} (t\sqrt{2})\right) \right )
\end{aligned}
```

Or as a MATLAB one-liner:

```matlab
@(x, n, s, T) 4/T * (x*sqrt(2)/2 + sum(exp(-s^2 * ((1:n).').^2/2 * (2*pi/T)^2) .* sin((1:n).' .* x * sqrt(2) * (2*pi/T)) ./ ((1:n).' * (2*pi/T))));
```

And comparing the two error functions, using ten sines in our approximation with a period of 10:

<p align="center">
    <img src="./figures/erf approx.png" width="75%"><br>
    <i>Figure 4: Comparing erf and Its Approximation</i>
</p>

But is this slower? Is it less accurate?  The one liner is most definitely always going to be slower, as it requires a whole bunch of memory be allocated and then worked through, so depending on the total number of sines you want to use, and the size of the input (e.g. if it's super large to test how fast it is), you can spend way too much time dealing with memory.  Or, you can loop through each sine instead; which, while not a convenient one-liner, seems to perform better.  Though, since the error function is typically approximated using a few polynomials as opposed to however each sine function is calculated, we can assume this approach will definitely be slower.

But how much?  Let's find out:

```matlab
>> x = linspace(-10, 10, 2^20);

>> tic; for i = 1:100, u = erf((x - 1) / 2); end, toc;
Elapsed time is 0.344611 seconds.

>> tic; for i = 1:100, u = baderf((x - 1) / 2, 10, 10); end, toc;
Elapsed time is 1.672149 seconds.

>> 1.672149 / 0.344611

ans =

    4.8523

>>  
```

Almost five times slower, that's not so bad!  We can, of course, speed that up by using fewer sines at the expense of accuracy, but we can never quite get as fast as the built in `erf`.  Oh well.  What about accuracy?  Well, the built in functions are pretty good with respect to that, but we can run some quick checks, say, a standard deviation out to see how various levels of approximation work out:

```matlab
>> for i = 1:10, err(i) = baderf(1, i, 10) - erf(1), end
>> err = table((1:10).', err.', 'VariableNames', ["nmax", "error"])

err =

  10x2 table

    nmax            error        
    ____    _____________________

      1        -0.154243549125533
      2       -0.0127841678455513
      3       0.00366666368192969
      4      0.000953672091459468
      5       7.1027478101171e-05
      6      1.62817545112937e-07
      7      -1.9881938284616e-07
      8      -8.0543699487734e-09
      9     -8.65575389141782e-11
     10      8.62643290133747e-13
```

Which, if you want to be able to quickly choose your accuracy, isn't so bad.  Assuming that the MATLAB built in error function is perfectly accurate at one standard deviation, and you only need six digits of accuracy, you can cut your run time down to only, like, three times slower if you choose six sines to use in your approximation instead of ten.  Of course, changing the period of the approximation also plays a big role.  After all, if you spread the function out, you'll probably have a worse accuracy for a given number of summed sines.  For example:

```matlab
>> for i = 1:10, err(i) = baderf(1, 10, i*4) - erf(1), end
>> err = table((4:4:40).', err.', 'VariableNames', ["T", "error"])

err =

  10x2 table

    T             error        
    __    _____________________

     4      0.00971564930514868
     8     4.52483606139253e-11
    12    -3.54251405987327e-09
    16     4.92961805265146e-07
    20     0.000186347314430124
    24      0.00131448958765756
    28      0.00294086686587258
    32      0.00269042300725186
    36     -0.00166862221701158
    40      -0.0109954993238769
```

Which is kind of fun.  Accuracy starts out really poor for a given number of summed sines at one standard deviation as adjacent periods of the function smoosh up against each other, but then as things stretch out again accuracy improves, but then quickly starts to deteriorate.  So yeah, we have multiple levers we can pull in order to adjust accuracy... yay...

<!-- ```math
\begin{aligned}
f(t - \mu) &= \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2\omega^2} \mathrm{e}^{-j\omega \mu} \mathrm{e}^{j\omega t} \; \mathrm{d}\omega \\

&= \frac{1}{2\pi} \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2\omega^2} \mathrm{e}^{j\omega (t-\mu)} \; \mathrm{d}\omega \\

&= \int_{-\infty}^\infty \mathrm{e}^{-\frac{1}{2} \sigma^2(2\pi f)^2} \cos(2\pi f (t-\mu)) \; \mathrm{d}f
\end{aligned}
``` -->
