import time
import sys

def onMultElemLine(A, B):
    n = len(A)
    C = [[0] * n for _ in range(n)]

    for i in range(n):  # Iterate over each row of A
        for k in range(n):  # For each element A[i][k]
            for j in range(n):  # Multiply by the row B[k]
                C[i][j] += A[i][k] * B[k][j]
    return C

def generateMatrix(n):
    # Generate a square matrix of size n x n filled with values
    return [[1 for _ in range(n)] for _ in range(n)]

def testOnMultLine(size):
    # Generate two matrices of size (size x size)
    A = generateMatrix(size)
    B = [[i + 1 for _ in range(size)] for i in range(size)]

    # Record the start time
    time1 = time.time()

    result = onMultElemLine(A, B)

    # Record the end time
    time2 = time.time()

    elapsed_time = time2 - time1

    # Display first 10 elements of the result matrix
    print("\nResult matrix:", file=sys.stderr)
    for i in range(1):
        for j in range(min(10, size)):
            print(result[i][j], end=" ", file=sys.stderr)
    print("", file=sys.stderr)

    print(f"{elapsed_time:.3f}")

    return elapsed_time

if __name__ == "__main__":
    size = int(sys.argv[1])
    testOnMultLine(size)