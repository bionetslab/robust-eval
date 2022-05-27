import time
from urllib.request import urlretrieve
import requests
import graph_tool as gt
import networkx as nx   # read_graphml from repotrial networks works with networkx 2.2 and not 2.5 which is the lastest one!!!


# get the network containing protein-protein with proper parameters via API
base_url = "https://api.nedrex.net"
submit_url = f"{base_url}/graph_builder"

data = {
    "nodes" : [],
    "edges" : ["protein_interacts_with_protein"],
    "concise" : True
}

print("Submitting request")
gbuild = requests.post(submit_url, json=data)
print(f"UID for job: {gbuild.json()}")
uid = gbuild.json()

while True:
    progress = requests.get(f"{base_url}/graph_details/{uid}")
    built = (progress.json()["status"] == "completed")
    if built:
        break
    print("Waiting for build to complete, sleeping for 10 seconds")
    time.sleep(10)

fname = "temp-PPI"
urlretrieve(f"{base_url}/graph_download_v2/{uid}/{fname}.graphml", f"{fname}.graphml")


G = nx.read_graphml(f"{fname}.graphml")
G_und = G.to_undirected()

node_list = set(G_und.nodes)
nodeAttr_list = {'geneName', 'taxid', 'type', 'displayName'}
for n in node_list:
    for attr in nodeAttr_list:
        if attr in G_und.nodes[n].keys():
            del G_und.nodes[n][attr]

edge_list = set(G_und.edges)
edgeAttr_list = {'memberOne', 'memberTwo', 'reversible', 'type', 'evidenceTypes'}
for e in edge_list:
    for attr in edgeAttr_list:
        if attr in G_und.edges[e].keys():
            del G_und.edges[e][attr]

network_name = "PPI-temp-graph.graphml"
nx.write_graphml(G_und, network_name)

gg = gt.load_graph(network_name)
gg.save("PPI.gt")