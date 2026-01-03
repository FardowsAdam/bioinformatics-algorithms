# Name: Fardows Adam
# ID: 4313204
# Email: 4313204@upm@edu.sa
# Assignment 4: Genome Assembly using De Bruijn Graph


def make_graph(reads, k):
    # Build the graph from reads
    graph = {}
    in_count = {}
    out_count = {}
    
    # Go through each read
    for r in reads:
        # Make all k-mers
        for i in range(len(r)-k+1):
            left = r[i:i+k-1]
            right = r[i+1:i+k]
            
            # Add to graph
            if left not in graph:
                graph[left] = []
            graph[left].append(right)
            
            # Count edges
            out_count[left] = out_count.get(left, 0) + 1
            in_count[right] = in_count.get(right, 0) + 1
    
    # Make sure all nodes are counted
    all_nodes = set(list(graph.keys()) + list(in_count.keys()))
    for node in all_nodes:
        if node not in out_count:
            out_count[node] = 0
        if node not in in_count:
            in_count[node] = 0
    
    return graph, in_count, out_count

def find_path(graph, in_count, out_count):
    # Find start node
    start = None
    for node in graph:
        if out_count[node] > in_count[node]:
            start = node
            break
    
    if start is None:
        start = list(graph.keys())[0] if graph else None
    
    if not start:
        return []
    
    # Make copy for walking
    temp_graph = {k: v[:] for k, v in graph.items()}
    
    # Walk through graph
    stack = [start]
    path = []
    current = start
    
    while True:
        # Check if we can go forward
        if current in temp_graph and temp_graph[current]:
            stack.append(current)
            # Go to next
            current = temp_graph[current].pop()
        else:
            # Add to path and go back
            path.append(current)
            if not stack:
                break
            current = stack.pop()
    
    # Reverse the path
    path.reverse()
    return path

def assemble(path):
    if not path:
        return ""
    
    # Build sequence
    result = path[0]
    for i in range(1, len(path)):
        result += path[i][-1]
    
    return result

# Main code
print("Genome assembly")
print("Enter k:")
k = int(input())
print("Enter reads (space separated):")
reads = input().split()
graph, in_cnt, out_cnt = make_graph(reads, k)
path = find_path(graph, in_cnt, out_cnt)
sequence = assemble(path)

print("Result:")
print(sequence)
