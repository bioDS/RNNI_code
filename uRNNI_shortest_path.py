__author__ = 'Lena Collienne'
#Computes exact shortest path in RNNI by using algorithm of Gavryushkin et al. (2018) for computing the RNNI graph and then computing distaces with BFS


from uRNNI_graph import *

tree_list = read_trees_from_file()

tree1 = string_to_list_of_sets(tree_list[0])
tree2 = string_to_list_of_sets(tree_list[1])


if is_tree(tree1) and is_tree(tree2):
    tree1 = convert_input_to_tree(tree1)
    tree2 = convert_input_to_tree(tree2)
    path = shortest_uRNNI_path(tree1, tree2)
    print ('A shortest path:')
    for tree in path[1]:
        tree = string_to_tree(tree)
        output = "["
        for i in tree:
            current_string = "{"
            for j in i:
                current_string = current_string + str(j) + ","
            current_string = current_string[:-1] + "}"
            output = output + current_string + ","
        output = output[:-1] + "]"
        print (output)
    print ('Distance:')
    print (path[0])
elif is_tree(tree2):
    print ('The first list does not represent a tree!')
else:
    print ('The second list does not represent a tree!')
