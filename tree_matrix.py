__author__ = '@gavruskin'
# This approach uses matrix representation of splits as in https://github.com/gavruskin/tauGeodesic
# For curvature computation, the set approach as implemented in tree.py is more relevant.


class Partition(object):
    # Partition is a matrix and a time coordinate.
    def __init__(self, matrix, coordinate):
        self.matrix = matrix
        self.coordinate = coordinate
        self.size = len(matrix)


class Tree(object):
    # Tree is a list of partitions.
    def __init__(self, partitions):
        self.partitions = partitions
        self.num_tips = partitions[0].size


def refines(partition1, partition2):
    # Returns 0 if the partitions are identical, 1 if partition 1 refines (cuts the tree closer to tips) partition 2,
    # -1 if partition2 refines partition1, and 2 if they are incompatible.
    # .coordinate are ignored!
    if partition1.size != partition2.size:
        return 2
    for i in partition1.matrix:
        if len(i) != partition1.size:
            return 2
    for i in partition2.matrix:
        if len(i) != partition1.size:
            return 2
    if partition1.matrix == partition2.matrix:
            return 0
    all_good = 0
    for i in range(partition1.size):
        all_good_tmp = 0
        for j in range(i, partition1.size):
            if not partition1.matrix[i][j] or partition2.matrix[i][j]:
                all_good_tmp += 1
        if all_good_tmp == partition1.size - i:
            all_good += 1
    if all_good == partition1.size:
        return 1
    all_bad = 0
    for i in range(partition1.size):
        all_bad_tmp = 0
        for j in range(i, partition1.size):
            if not partition2.matrix[i][j] or partition1.matrix[i][j]:
                all_bad_tmp += 1
        if all_bad_tmp == partition1.size - i:
            all_bad += 1
    if all_bad == partition1.size:
        return -1
    else:
        return 2


def refines_time(partition1, partition2):
    # Same as refines but takes times into account
    if refines(partition1, partition2) == 0 and partition1.coordinate == partition2.coordinate:
        return 0
    elif refines(partition1, partition2) == 1 and partition1.coordinate < partition2.coordinate:
        return 1
    elif refines(partition1, partition2) == -1 and partition1.coordinate > partition2.coordinate:
        return -1
    else:
        return 2


def are_compatible(tree, partition):
    # Returns 1 if partition is compatible with tree, otherwise 0.
    for tree_part in tree.partitions:
        if refines(tree_part,partition) == 2:
            return 0
    return 1


# TODO: For future, implement are_compatible_time(tree, partition)


def sort(tree):
    # Sorts partitions in tree under the refinement relation.
    # Cold be avoided to speed things up depending on the purpose.
    # Implements bubble sort.
    done = False
    while not done:
        done = True
        for i in range(len(tree.partitions)-1):
            if refines(tree.partitions[i], tree.partitions[i+1]) == 2:
                print("The input has to be a tree. Your input has incompatible partitions.")
                return
            elif refines(tree.partitions[i], tree.partitions[i+1]) == 0:
                print("The input tree has identical partitions.")
                return
            elif refines(tree.partitions[i], tree.partitions[i+1]) == -1:
                tree.partitions[i], tree.partitions[i+1] = tree.partitions[i+1], tree.partitions[i]
                done = False


def add_partition(tree, partition):
    if not are_compatible(tree, partition):
        print("You are trying to add a partition which is not compatible with the tree.")
        return
    for part in tree.partitions:
        if refines(part, partition) == 0:
            print("You are trying to add a partition which is already present in the tree.")
            return
    else:
        tree.partitions.append(partition)
        sort(tree)
