# For R Users

Here is a gentle introduction to Julia and `GLV` destined to R users.

## Download Julia

To Download and install Julia, please go to the official [Julia website](https://julialang.org/downloads/).

Once you have installed Julia, you are ready to go.
The best to follow this small guide is to open a Julia terminal 
and type the commands as you read.
This way you will get a better understanding of how Julia works.

## Julia Makes You A Better Programmer

At first, Julia can be a *very* intimidating language.
It is often more strict than R and may throw many errors, especially in the beginning.
But don’t worry, it will get easier with time.
In fact, when Julia throws errors, it’s often for your own good.

For example, in R, it is valid to add two vectors of different lengths like this

```r
c(1, 2, 3, 4) + c(1, 2)
```

However, the equivalent code in Julia will throw an error

```julia
[1, 2, 3, 4] + [1, 2]
```
::: tip

In `R` vectors are written using `c()` and in `Julia` they are written using the brackets `[]`.
For an exhaustive list of noteworthy differences between `R` and `Julia`, see [this page](https://docs.julialang.org/en/v1/manual/noteworthy-differences/#Noteworthy-differences-from-R).

:::

The error occurs because the operation **ambiguous**.
When the second vector runs out of elements, it’s unclear whether it should repeat the smaller vector or simply stop adding. 
For example, should the result be `[2, 4, 4, 6]` (by repeating the smaller vector), 
or `[2, 4, 3, 4]` (by adding nothing when the smaller vector ends)?

Since this behaviour is non-trivial, Julia asks you to be explicit. 
If you want the second vector to repeat, say it

```julia
[1, 2, 3, 4] + repeat([1, 2], 2)
```

This not only ensures that your code behaves as intended, 
but it also makes it clearer and easier for others (and yourself) to understand.

This is was a simple example explaining that if Julia is punishing you with errors,
it is often to make you write better code.
And the other main reason is to make your code run faster.



