
# coding: utf-8

# In[26]:

import networkx as nx



inputFile = "/Users/soheil.danesh/Documents/PROJECTS/Mediator Complex Functional Differences in Homologs/BIOGRID-ORGANISM-3.4.143.tab/BIOGRID-ORGANISM-Homo_sapiens-3.4.143.tab"
rdd = sc.textFile(inputFile)
rdd.take(5)

proteinRelationsTable = rdd.collect()

        
        
    
    

rdd = rdd.map(lambda row : row.split('\t') )

#rdd.coalesce(1).saveAsTextFile('/Users/soheil.danesh/Documents/PROJECTS/hubClique/biogrid-homoSapiens-ppi.csv')

proteinRelationsTable = rdd.collect()

for i in range(0, 10):
    print(proteinRelationsTable[i])

#read input file



# In[27]:

#make graph
#a node consists of an ngram(term) and its weight(Statistical)

def makeGraph(relationsTable, nodeColumns=[2,3]): #todo maybe later add a weight column

    if(len(nodeColumns) != 2):
        raise 'nodeColumns is supposed to have two members each being the index of the column in the table where the nodes of the graph are'

    G = nx.Graph()
    for row in relationsTable:
        node1 = row[nodeColumns[0]]
        node2 = row[nodeColumns[1]]
        G.add_edges_from([(node1,node2)])  
        
    return G
            

#for source med12

#find neightbors

#for each neighbor that has no weight assignments yet assign a new neigh


# In[28]:

graph = makeGraph(proteinRelationsTable)


# In[29]:

from collections import deque

def updateNodeWeightsMap( node, weight, nodeWeightsMap):
    firstTimeSeeingNode = True 
    if(node in nodeWeightsMap):
        nodeWeightsMap[node] += weight
        firstTimeSeeingNode = False
    else:
        nodeWeightsMap[node] = weight
        
    return firstTimeSeeingNode


#the way blast works is as the activation diffuses, we make sure each node emenates only once by adding only weights that are seeing for the first time to the sources queue and dequeueing after one emanation
#btw difusion happens by dividing actvation reached a node by the node's degree and sending it to its neightbor, with an optional decay factor also multiplied into each weight, 
def spreadingActivationSourceWeight(graph, sourceNodeName, weight=1.0, decayFactor=1.0):
    
    nodeWeightsMap = {}
    
    sourcesq = deque([ (sourceNodeName, weight) ])
    
    while(len(sourcesq)>0):
        nodeWeight = sourcesq.popleft()
        nodeName = nodeWeight[0]
        
        if(graph.degree(nodeName) == 0):
            raise 'node has degree of 0, develop this code branch'
            
        emenatedWeight = nodeWeight[1] * decayFactor / float(graph.degree(nodeName))
        for neighbor in graph[nodeName]:
            #add the emenated weight to the neighbor and also see whether it is the first time we are seeing this neighbor or not
            firstTimeSeeingNode = updateNodeWeightsMap(neighbor, emenatedWeight, nodeWeightsMap)
            if(firstTimeSeeingNode):# if it's the first time we are seeing this node make sure it emenates in the next iteration
                sourcesq.append( (neighbor, emenatedWeight) )
                
    #put node-weight pairs in nodeWeightMap into a sorted array as well
    nodeWeightsArr = []
    for node in nodeWeightsMap:
        nodeWeightsArr.append( (node, nodeWeightsMap[node]) )
        
    nodeWeightsArr = sorted(nodeWeightsArr, key=lambda nodeWeight: nodeWeight[1])
        
    return (nodeWeightsArr, nodeWeightsMap)


    #TODO?
    #def normalizeNodeWeightArr():
        
    
    
    


# In[30]:

(nodeWeightsArr, nodeWeightsMap) = spreadingActivationSourceWeight(graph, 'MED12')
for nodeWeight in nodeWeightsArr:
    print(nodeWeight)


# In[55]:

(nodeWeightsArr, nodeWeightsMap) = spreadingActivationSourceWeight(graph, 'MED12L')
def printNodeWeightsArr(nodeWeightsArr, limit=-1):
    if(limit == -1):
        limit = len(nodeWeightsArr)
        
    rank = 0
    for nodeWeight in nodeWeightsArr[0:limit]:
        rank += 1
        print("rank = "+ str(rank) + " nodeWeight = "+ str(nodeWeight))
printNodeWeightsArr(nodeWeightsArr,20)


# In[63]:



#vector normalization, addition, subtraciton and dot product using spark

#vector is a rdd of tuples like (u'CDK8', 0.13727614584604667) ie [(u'CDK8', 0.13727614584604667), .....]
    #bookmark: make these functions, then for thyroid cancer:
    #after normalizing vectors
    # add the vectors of the source nodes to find a list of likely candidate genese that may also be involved
    
    #for med complex subtract the vecotre of individual homolougs e.g. med12 - med12L

    

def normalizNodeWeightsArr(nodeWeightsArr):
    nodeWeightsRdd = sc.parallelize(nodeWeightsArr)
    sumNormalizedRdd = normalizNodeWeightsRdd(nodeWeightsRdd)
    sumNormalizedSortedNodeWeightsArray = sumNormalizedRdd.map(lambda (a,b) : (b,a)).sortByKey(ascending=False).map(lambda (a,b) : (b,a)).collect()
    return sumNormalizedSortedNodeWeightsArray

def normalizNodeWeightsRdd(nodeWeightsRdd):
    weightsSum = nodeWeightsRdd.map(lambda nodeWeight: nodeWeight[1]).sum()
    sumNormalizedRdd = nodeWeightsRdd.map(lambda nodeWeight: (nodeWeight[0], float(nodeWeight[1]) / float(weightsSum) )) 
    return sumNormalizedRdd

def normalizeAndAddNodeWeightArrs(nodeWeightsArrs):
    return normalizeAndOperateOnNodeWeightArrs(nodeWeightsArrs, operation='add')

def normalizeAndSubtractNodeWeightArrs(nodeWeightsArrs):
    return normalizeAndOperateOnNodeWeightArrs(nodeWeightsArrs, operation='subtract')


def normalizeAndOperateOnNodeWeightArrs(nodeWeightsArrs, operation='add'):
    resultsRdd = normalizNodeWeightsRdd(sc.parallelize(nodeWeightsArrs[0]))
    for arr in nodeWeightsArrs:
        resultsRdd.union(normalizNodeWeightsRdd(sc.parallelize(arr)))
        if(operation == 'add'):
            resultsRdd = resultsRdd.reduceByKey(lambda a,b: a+b)
        elif( operation == 'subtract'):
            resultsRdd = resultsRdd.reduceByKey(lambda a,b: a-b)
            
    
    sortedresultsArr = resultsRdd.map(lambda (a,b) : (b,a)).sortByKey(ascending=False).map(lambda (a,b) : (b,a)).collect() #from http://stackoverflow.com/questions/24656696/spark-get-collection-sorted-by-value/40708098#40708098 
    return sortedresultsArr
        
        
    
#def subtractVectors(): tood bookmark

    
    


# In[64]:

(nodeWeightsArr_med12, nodeWeightsMap_med12) = spreadingActivationSourceWeight(graph, 'MED12')
(nodeWeightsArr_med12L, nodeWeightsMap_med12L) = spreadingActivationSourceWeight(graph, 'MED12L')

med12LMinusMed12 = normalizeAndSubtractNodeWeightArrs([nodeWeightsArr_med12L, nodeWeightsArr_med12])

printNodeWeightsArr(med12LMinusMed12,20)


# In[66]:

med12MinusMed12L = normalizeAndSubtractNodeWeightArrs([nodeWeightsArr_med12, nodeWeightsArr_med12L])
printNodeWeightsArr(med12MinusMed12L,50)


# In[ ]:




# In[59]:

#sumNormalizedSortedNodeWeightsArray = normalizNodeWeightsArr(nodeWeightsArr)
#print(sumNormalizedSortedNodeWeightsArray)


# #thyroid cancer genes
# Text of email from Vala about genes involved in thyroid cancer,
# eblow I've studied BRAF and TERT and show that the TERT promoters are all ranked in the top 50 or so ranked genes. also RAF1 seems involved in thyroid cancer in a 2015 paper:
# miR-195 is a key regulator of Raf1 in thyroid cancer
# (https://www.ncbi.nlm.nih.gov/pubmed/26527888)
# 
# BRAF mutation co- occurring with a TERT promoter, PIK3CA, TP53, or AKT1 mutation (155–157,577). Such a combination of several mutations is seen in a much smaller fraction of PTC as compared with a 40%–45% incidence of BRAF mutations and is expected to serve as a more specific marker of unfa- vorable outcomes of PTC.
# 
# The combination of a TERT mutation and a BRAF mutation within the same tumor was associated with a high risk of structural disease recurrence 
# 
# seven-gene mutational panel (BRAF, RAS, PAX8/PPARc, RET/PTC) would be more likely to demonstrate similar performance to that of a nonpregnant population.
# 
# Mutational testing for BRAF in AUS/FLUS samples has high specificity for cancer, but low sensitivity (198,199). Testing for a panel of mutations (BRAF, NRAS, HRAS, KRAS, RET/PTC1, RET/PTC3, PAX8/PPARc) offers a significantly higher sensitivity of 63%–80% ( 
# 

# In[ ]:

#thyroid cancer genes

(braf_nodeWeightsArr, braf_nodeWeightsMap) = spreadingActivationSourceWeight(graph, 'BRAF')
#print(braf_nodeWeightsArr)
(tert_nodeWeightsArr, tert_nodeWeightsMap) = spreadingActivationSourceWeight(graph, 'TERT')

braf_and_tert_nodeWeightArr = normalizeAndAddNodeWeightArrs([braf_nodeWeightsArr, tert_nodeWeightsArr])

    


# In[56]:

printNodeWeightsArr(braf_and_tert_nodeWeightArr, 50)


# In[ ]:




# In[ ]:


#This code was done while experimenting to make the method 
#simple diffusion
#rules: each source starts with a labled weight score of 1. The lable could be the protin itself (To find related proteins) or a known funciton of it (to assign functions to proteins)
from collections import deque


sourcesq = deque([ ('MED12', 1) ])
nodeWeightsMap = {}
decayFactor = 0.95

def addToNodeWeightsMap( node, weight ):
    firstTimeSeeingNode = True 
    if(nodeWeight[0] in nodeWeightsMap):
        nodeWeightsMap[nodeWeight[0]] += nodeWeight[1]
        firstTimeSeeingNode = False
    else:
        nodeWeightsMap[nodeWeight[0]] = nodeWeight[1]
        
    return firstTimeSeeingNode


while(len(sourcesq)>0):
    nodeWeight = sourcesq.popleft()
    nodeName = nodeWeight[0] 
    nodesWeight = nodeWeight[1] 
    for neighbor in graph[nodeName]:
        firstTimeSeeingNode = addToNodeWeightsMap(nodeName, nodesWeight*decayFactor)
        if(firstTimeSeeingNode):
            sourcesq.append(nodeName, nodesWeight*decayFactor)
        
    



sources = [('MED12', 1)]
decayFactor = 0.95

nodeWeights = {} #keeps track of the list of nodes that have already emenated their activaiton
alreadyEmenated = [] #holds nodes that have already emenated their weights, a node can only emenate weights once for each weight label, to avoid flowback and oscilations



for nodeWeight in sourcesq:
    #alreadyEmenated.append(node)
    for neighbor in graph[nodeWeight[0]]:
        #if(neighbor not in alreadyEmenated)
            if(neighbor in nodeWeights):
                nodeWeights[neighbor] += nodeWeight[0] * decayFactor 
            else:
                nodeWeights[neighbor] = nodeWeight[0] * decayFactor 
        
    


# In[12]:

print(graph['MED12'])


# In[20]:

import community
parts = community.best_partition(graph)
#print(part)
#print("-------")

#parts is an array where keys are nodes, values are paritions
partitions = {}
for node in parts: 
    partition = parts[node]
    if( partition in partitions):
        partitions[partition].append(node)
    else:
        partitions[partition] = [node]
        
for key in partitions:
    print("-------------")
    print(partitions[key])
    print("-------------")
        


# In[21]:

graph['MED12']


# In[22]:

graph['MED12L']


# In[ ]:



