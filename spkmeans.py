import math
import sys
import numpy as np


def get_CMD_input():
    if len(sys.argv) == 4:
        temp_k = sys.argv[1]
        goal = sys.argv[2]
        file_name = sys.argv[3]
        if (sys.argv[1]).isdigit():
            k = int(temp_k)
            if k < 1:
                print("Invalid input!")
                sys.exit(1)
    elif len(sys.argv) == 3:
        goal = sys.argv[1]
        file_name = sys.argv[3]
        k = 0
    if (len(sys.argv) < 3) or (len(sys.argv) > 4):
        print("Invalid Input!")
        sys.exit(1)
    return goal, k, file_name


goal, k, file_name = get_CMD_input()
