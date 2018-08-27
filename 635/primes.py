
class Cache:
    invfac, invfac2 = 0, 0

    def __init__(self):
        pass

def primes_dict_from_file(maximal):
    primes = dict()
    with open('primes/total.txt', 'r') as filehandle:
        for line in filehandle:
            for s in line.split('\t'):
                if(int(s) > maximal):
                    return primes
                primes[int(s)] = Cache()

        return primes

# print(len(primes_dict_from_file(10000000000)))
