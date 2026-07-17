using Random
using Statistics

function testMethod(array::Vector{Float64})
    arr = deepcopy(array)
    t = 0.0
    r = 0.0
    s = 0.0
    c = 0.0
    i = 1
    del = []
    while length(arr) > 0
        del = []
        s = sum(arr)
        t += randexp() / s
        # println((t,s, length(arr)))

        for _ in 1:length(arr)
            r = rand() * s
            # println((t, r, s, length(arr)))
            i = 1
            c = arr[i]
            while c < r
                i += 1
                c += arr[i]
            end
            push!(del, i)
        end

        del = sort(unique(del))

        for j in length(del):-1:1
            deleteat!(arr, del[j])
        end
    end

    return t
end

function refMethod(arr::Vector{Float64})
    return randexp() / minimum(arr)
end

function test(l,iters)
    arr = rand(l)*0.1
    ref = mean([refMethod(arr) for _ in 1:iters])
    test = mean([testMethod(arr) for _ in 1:iters])
    # println("Array length: ", length(arr))
    # println("Array minimum: ", minimum(arr))
    # println("Iterations: ", iters)
    # println("Reference: ", ref)
    # println("Test: ", test)
    println(length(arr), " items; ",ref,", ",test,"; Percent error: ", 100 * (test - ref) / ref, "%")
end

function testArr(arr,iters)
    ref = mean([refMethod(arr) for _ in 1:iters])
    test = mean([testMethod(arr) for _ in 1:iters])
    println(length(arr), " items; ",ref,", ",test,"; Percent error: ", 100 * (test - ref) / ref, "%")
end

# test(ones(1000), 10000)
# test(rand(2)*2, 10000)
# for i in 1:10
#     test(i, 1000000)
# end
#
for i in 1:20
    testArr([1.,i], 1000000)
end
