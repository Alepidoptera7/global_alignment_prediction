
"""The purpose of this algorithm is to provide all prime numbers up to a user provided limit, n."""

#initial conditions
n = 77
odd_list = [2]

for i in range(2, n):
    #generates all odd numbers, cutting number of computations in half
    if i % 2 != 0:
        odd_list.append(i)

prime_list = []
for i in odd_list:
    division_check_list = []
    for j in range(1, n):
        if i % j == 0:
            #i is divisible by this number, j
            division_check_list.append(i)

    if len(division_check_list) == 2:
        #if the number of integers that a number is divisible by is exactly one and itself,
        # the length of the division check list will be exactly 2
        #prime condition for i
        prime_list.append(i)

print("prime_list", prime_list)