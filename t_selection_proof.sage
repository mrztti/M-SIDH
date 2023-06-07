


def prove_t(sec_pam,condition, results, restart_with_t=None):

    # sec_pam is lambda security parameter
    # we start at t=2lambda
    if restart_with_t == None:
        t = 2*sec_pam
    else:
        t = restart_with_t
    
    # sample the 2t smallest primes
    sample_length = 2*t
    primes = Primes()
    primes_list = [primes.unrank(0) ** 2]
    # collect the primes
    for i in range(1, sample_length):
        primes_list.append(primes.unrank(i))
    
    primes_A = primes_list[::2]
    primes_B = primes_list[1::2]
    
    assert len(primes_A) == len(primes_B)
    assert len(primes_A) == t

    B = prod(primes_B)

    n = 0  
    while B <= prod(primes_A[(n):]) ** 2:
        n+=1
        # because we add one after the loop condition, we add an extra index
        # this is correct, as the index n in the paper starts at 1 not 0

    res = results + [(t, n)]
    # condition is a function (Int, Int) -> Bool
    if condition(sec_pam, (t - n + 1)):
        
        return prove_t(sec_pam, condition, res, restart_with_t=t+1)
    
    return t, n , res


if __name__ == "__main__":

    paper_condition = lambda a,b : a < b
    correct_condition = lambda a,b : a > b

    print("Condition from paper: lambda < t-n+1")
    print("Correct condition: lambda > t-n+1")


    t1, _, _  = prove_t(128, paper_condition, [])
    t2, _, _  = prove_t(128, correct_condition, [])

    # correct rank is 2t + 1 since rank in python start at 0 and we take sample length 2t
    print(f"For lambda = 128, the paper condition yields a max_prime of rank {2*t1 + 1}")
    print(f"For lambda = 128, the correct condition yields a max_prime of rank {2*t2 + 1}")
    expected = 571
    print(f"Expected rank: {expected}")


    t1, _, _  = prove_t(192, paper_condition, [])
    t2, _, _  = prove_t(192, correct_condition, [])

    print(f"For lambda = 192, the paper condition yields a max_prime of rank {2*t1 + 1}")
    print(f"For lambda = 192, the correct condition yields a max_prime of rank {2*t2 + 1}")
    expected = 851
    print(f"Expected rank: {expected}")


    t1, _, _  = prove_t(256, paper_condition, [])
    t2, _, _  = prove_t(256, correct_condition, [])

    print(f"For lambda = 256, the paper condition yields a max_prime of rank {2*t1 + 1}")
    print(f"For lambda = 256, the correct condition yields a max_prime of rank {2*t2 + 1}")
    expected = 1131
    print(f"Expected rank: {expected}")


