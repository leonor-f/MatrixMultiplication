import time
import sys

def onMult(A, B):
    # Matrix multiplication: C = A * B
    n = len(A)  # Assuming square matrices (n x n)
    m = len(B[0])  # Number of columns in B
    result = [[0] * m for _ in range(n)]  # Initialize the result matrix

    for i in range(n):
        for j in range(m):
            temp = 0
            for k in range(n):
                temp += A[i][k] * B[k][j]
            result[i][j] = temp
    return result

def generateMatrix(n):
    # Generate a square matrix of size n x n filled with values
    return [[1 for _ in range(n)] for _ in range(n)]

def testOnMult(size):
    # Generate two matrices of size (size x size)
    A = generateMatrix(size)
    B = [[i + 1 for _ in range(size)] for i in range(size)]

    # Record the start time
    time1 = time.time()

    result = onMult(A, B)

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
    testOnMult(size)