from copy import deepcopy
from typing import Union


def solve(arr_a, arr_b: Union[list, tuple]) -> tuple:
    """Find the solution of the SLAE.
    :param arr_a: The coefficient matrix of equations
    :param arr_b: The column vector of right-hand-sides of equations
    :return: X, the solution vector
    """
    det = get_det_if_valid(arr_a, arr_b)

    arr_x = tuple(
        get_det(arr_replace_column(arr_a, arr_b, i)) / det
        for i in range(len(arr_a))
    )

    return arr_x


def get_det_if_valid(arr_a, arr_b: Union[list, tuple]) -> Union[int, float]:
    """Get the determinant of an array A, if arrays are correct.
    :param arr_a: The coefficient matrix of equations
    :param arr_b: The column vector of right-hand-sides of equations
    :return: The determinant of an array A
    """
    if not (arr_a and arr_b):
        raise ValueError(
            "The matrix or the column vector of right-hand-sides of the "
            "equations is empty!"
        )

    size = len(arr_a)
    if len(arr_b) != size:
        raise ValueError(
            "The size of the column vector of right-hand-sides of the "
            "equations must be equal to the size of the matrix!"
        )

    if size != 1:
        for arr in arr_a:
            if not isinstance(arr, list) or len(arr) != size:
                raise ValueError("The matrix must be quadratic!")

    det = get_det(arr_a)
    if not det:
        raise ValueError("The matrix determinant must not be 0!")

    return det


def get_det(arr) -> Union[int, float]:
    """Get the determinant of an array.
    :param arr: An array
    :return: The array determinant
    """
    if len(arr) == 1:
        return arr[0]

    det = 0
    for i, j in enumerate(arr[0]):
        det += j * ((-1) ** i) * get_det(arr_remove(arr, i))

    return det


def arr_remove(arr, index: int):
    """Remove the 1st row and the index-th column from an array.
    :param arr: An array
    :param index: An index
    :return: A changed copy of the array
    """
    res = deepcopy(arr[1:])

    for i in res:
        i.pop(index)

    if len(res) == 1:
        return res[0]
    return res


def arr_replace_column(arr, col: Union[list, tuple], index: int):
    """Replace the index-th column from an array by a column vector.
    :param arr: An array
    :param col: A column vector
    :param index: An index
    :return: A changed copy of the array
    """
    if len(arr) == 1:
        return col

    res = deepcopy(arr)

    for i, j in enumerate(col):
        res[i][index] = j

    return res

array_a = [
  [2, 5, 4],
  [1, 3, 2],
  [2, 10, 9],
]
array_b = [30, 150, 110]

print(solve(array_a, array_b))

