__author__ = "Lena Collienne"

#Computes caterpillar path and distance between two caterpillar trees that need to be given in a txt files in cluster representation [{},{},..] (two rows, each one caterpillar tree)

from caterpillar import *

def main():
    tree_list = read_trees_from_file()

    tree1 = string_to_list_of_sets(tree_list[0])
    tree2 = string_to_list_of_sets(tree_list[1])

    if is_caterpillar(tree1) and is_caterpillar(tree2):
        path = shortest_caterpillar_path(tree1,tree2)
        print("A shortest caterpillar path:")
        for i in path:
            print(i)
        print("Distance between the given trees:")
        print(len(path)-1)
    elif is_caterpillar(tree1) == False and is_caterpillar(tree2) == False:
        print("Neither of both input trees is a caterpillar tree")
    elif is_caterpillar(tree1):
        print ("The second tree is no caterpillar tree")
    elif is_caterpillar(tree2):
        print ("The first tree is no caterpillar tree")

if __name__ == "__main__":
    main()
