# Flexle

Fast, dynamically weighted random sampling.

## API

| Function/Operation | Description |
|---|---|
| `FlexleSampler(weights)` | Create a `FlexleSampler` object from a collection of `weights`. | 
| `sampler[i]`[^1] | Get the weight of the element at index `i` in `sampler`. |
| `sampler[i] = w`[^2] | Update the weight of the element at index `i` in `sampler` to be equal to `w`. |
| `push!(sampler, w)` | Add a new element of weight `w` to `sampler`, placing it at the end of `sampler.weights`. |
| `deleteat!(sampler, i)` | Remove element `i` from `sampler`, shifting the index of every subsequent element over to fill the gap. |
| `sample(sampler)` | Take a single random sample from `sampler`, returning the index of the element sampled. |


[^1]: Implemented as `getindex(sampler, i)`.
[^2]: Implemented as `setindex(sampler, w, i)`. When `setindex` is called explicitly, it also returns the difference between the new and old weight of element `i`.
