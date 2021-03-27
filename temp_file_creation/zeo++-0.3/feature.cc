//#include "network.h"
#include <cmath>

#include "feature.h"
#include "general.h"

using namespace std;

int FEATURE::createFeatures(ATOM_NETWORK *uc,VORONOI_NETWORK *origVNet, DIJKSTRA_NETWORK *origDNet, std::fstream &output, int initIndex,char *name){
    nfeatures = 0;
    features = std::vector<int> ();
    feature_node_mapping = std::vector<int> ();
    segmentChannel(uc);
    if(nfeatures == 0){
        std::cerr << "Error occurred during channel segmentation. No features were found." << "\n"
        << "Exiting..." << "\n";
        exit(1);
    }
    
    std::vector< std::vector<int> > nodesForFeatures = std::vector < std::vector<int> > (nfeatures, std::vector<int> ());
    std::vector<FEATURE> featChannels = std::vector<FEATURE> ();
    for(unsigned int i = 0; i < feature_node_mapping.size(); i++){
        nodesForFeatures[feature_node_mapping[i]].push_back(reverseIDMappings.find(i)->second);
    }
    
    
    int zeroVecs [3][3] = {{0, 0, 0},{ 0, 0, 0}, {0, 0, 0}};
    
    for(int i = 0; i < nfeatures; i++){
        DIJKSTRA_NETWORK featNet;
        DIJKSTRA_NETWORK::filterDnetEdges(nodesForFeatures[i], origVNet, &featNet);
        std::vector<bool> infoStorage; std::vector<CHANNEL> channels;
        findChannels(&featNet, &infoStorage, &channels);
        int dim;
        if(channels.size() != 0){
            dim = channels[0].dimensionality;
            featChannels.push_back(FEATURE(nodesForFeatures[i], origDNet, dim, channels[0].basis));
        }
        else{
            dim = 0;
            featChannels.push_back(FEATURE(nodesForFeatures[i], origDNet, 0, zeroVecs));
        }
        featChannels[i].writeToVMD("feature", initIndex + i, output);
        std::cout << "@@ " << name << " Feature " << initIndex + i << " volume = " << featChannels[i].calculateVolume() << "\n";
        featChannels[i].calculateSurfaceArea(uc,&nodes,&feature_node_mapping,i,name,initIndex + i);
    }
    return nfeatures;
}

int FEATURE::createSegments(ATOM_NETWORK *uc,VORONOI_NETWORK *origVNet, DIJKSTRA_NETWORK *origDNet, std::fstream &output, int initIndex){
    nsegments = 0;
    segments = std::vector<int> ();
    segmentChannel(uc);
    if(nsegments == 0){
        std::cerr << "Error occurred during channel segmentation. No segments were found." << "\n"
        << "Exiting..." << "\n";
        exit(1);
    }
    
    std::vector< std::vector<int> > nodesForSegments = std::vector < std::vector<int> > (nsegments, std::vector<int> ());
    std::vector<CHANNEL> segChannels = std::vector<CHANNEL> ();
    for(unsigned int i = 0; i < segments.size(); i++){
        nodesForSegments[segments[i]].push_back(reverseIDMappings.find(i)->second);
    }
    
    int zeroVecs [3][3] = {{0, 0, 0},{ 0, 0, 0}, {0, 0, 0}};
    
    for(int i = 0; i < nsegments; i++){
        DIJKSTRA_NETWORK segNet;
        DIJKSTRA_NETWORK::filterDnetEdges(nodesForSegments[i], origVNet, &segNet);
        std::vector<bool> infoStorage; std::vector<CHANNEL> channels;
        findChannels(&segNet, &infoStorage, &channels);
        int dim;
        if(channels.size() != 0){
            dim = channels[0].dimensionality;
            segChannels.push_back(CHANNEL(nodesForSegments[i], origDNet, dim, channels[0].basis));
        }
        else{
            dim = 0;
            segChannels.push_back(CHANNEL(nodesForSegments[i], origDNet, 0, zeroVecs));
        }
        segChannels[i].writeToVMD("segment", initIndex + i, output);
    }
    return nsegments;
}

void FEATURE::segmentChannel(ATOM_NETWORK *uc){
int i;

cout << "Current channel has " << nodes.size() << " nodes" << "\n";

for(i=0;i<(int)nodes.size();i++)
   {
   segments.push_back(-1);
   };
nsegments=0;

while(segment_findMaxNode()>-1)
   {
//   segment_growSegment_start(uc,segment_findMaxNode());
   segment_distBasedSegmentation(uc);
   };

cout << "Initial number of segments (minima) " << nsegments << "\n";

/*
fstream output1;
output1.open("segments.zvis",fstream::out);
segment_saveVis(output1);
output1.close();
*/

//segment_identify_connections(); // analyzes connectivity of each segment to find segment's 'faces' - connection between segments 



//merge_segments(); // analyze segments and merge if neccesary

//merge_nothing(); // this function copies segments into features

merge_singlenode_segments();


//calculateVolumeForAllFeatures(); // calculate volume of each segment


/* Print results */
cout << "\n" << "Segment info (ID - r):"<<"\n";
for(i=0;i<(int)segment_max_r.size();i++)
   {
   cout << i<< "   " << segment_max_r[i]  << "\n";
   };
cout << "\n" << "Segment connection info (from to radii merge_stat):"<<"\n";

/*
for(i=0;i<(int)segment_conn.size();i++)
   {
   cout << segment_conn.at(i).from_seg << "   " << segment_conn.at(i).to_seg << "   " <<segment_conn.at(i).max_radius << "   " << segment_conn.at(i).merged << "\n";
   };

cout << "\n" << "Total number of segment connecting edges :  "<< segment_conn_fulllist.size()  <<"\n";

cout << "\n" << "Segment info vs features (segment ID - feature ID):"<<"\n";
for(i=0;i<nsegments;i++)
   {
   cout << i<< "   " << segment_merge_list[i]  << "\n";
   };
*/
cout << "\n" << "Features info (ID - r - volume):"<<"\n";
for(i=0;i<(int)feature_max_r.size();i++)
   {
//   cout << i<< "   " << feature_max_r[i]  << "   " << feature_volume[i] << "\n";
   cout << i<< "   " << feature_max_r[i]   << "\n";
   };
/*
cout << "\n" << "Features connection info (from to radii ):"<<"\n";
for(i=0;i<(int)feature_conn.size();i++)
   {
   cout << feature_conn.at(i).from_seg << "   " << feature_conn.at(i).to_seg << "   " <<feature_conn.at(i).max_radius << "\n";
   };
*/

} // end of segmentChannel

/* find node with largest max_radius that is not assigned to any segment yet */
int FEATURE::segment_findMaxNode(){
int i;
int nd=-1;
double min=0;
for(i=0;i<(int)nodes.size();i++)
   {
   if(segments[i]==-1&&nodes.at(i).max_radius>min) 
       {
       nd=i;
       min=nodes.at(i).max_radius;
       };
   };
return nd;
} // end of segment_findMaxNode()



/* NEW segmentation rountine based on distance criterion */
/* this function needs to have a pointer to unit cell to run distance calcluations
 * please call setUnitCellPointer() prior to this function */
void FEATURE::segment_distBasedSegmentation(ATOM_NETWORK *uc)
{
UnitCellPtr=uc;

int cur_maxId; // id of current max node
double cur_r;
int i;

double scale_factor=1.1;

if(UnitCellPtr==NULL) {cout<< "no distance pointer set";abort();};


do{
/* clear variable */
segments.clear();
segment_max_r.clear();
for(i=0;i<(int)nodes.size();i++)
   {
   segments.push_back(-1);
   };
nsegments=0;

/* adjust scale factor */
scale_factor-=0.05;

if(scale_factor==0.5)
    {
    cerr << "Segmentation failed. Aborting" << "\n";
    abort();
    };

/* perform segmentation */
while(segment_findMaxNode()>-1)
  {
  cur_maxId=segment_findMaxNode();
  cur_r=nodes.at(cur_maxId).max_radius;
  segments[cur_maxId]=nsegments; // assign new segment number to a nw segment
  segment_max_r.push_back(nodes.at(cur_maxId).max_radius); //store maximum r in a segment

  // now assign all nodes with close distance to the current node
  for(i=0;i<(int)nodes.size();i++)
     {
     if(UnitCellPtr->calcDistanceXYZ(nodes.at(cur_maxId).x,nodes.at(cur_maxId).y,nodes.at(cur_maxId).z,nodes.at(i).x,nodes.at(i).y,nodes.at(i).z)<scale_factor*cur_r)
        {
        segments[i]=nsegments;
        };
     };

  nsegments++;
  };

}while(validateSegmentation()==0);
//}while(scale_factor>0.5);

cout << " no segments after segment_distBasedSegmentation: " << nsegments << "\n"; 

} // end of segment_distBasedSegmentation()


/* this function validates current segmentation
 * valid segments are build from interconnected nodes only */

int FEATURE::validateSegmentation(){
int i;
int valid=1;
for(i=0;i<nsegments;i++)
   {
   if(validateSegment(i)==0)
     {
     valid=0;
     break;
     };

   }

return valid;

} // end of CHANNEL::validateSegmentation()



int FEATURE::validateSegment(int segid){
int i;
int id,conn_to;
int nvisited=0;
int valid=1;
vector <int> nodelist;
vector <int> visited;

vector <int> heap;

/* count nodes and init heap */
for(i=0;i<(int)nodes.size();i++)
   {
   if(segments[i]==segid) 
      {
      nodelist.push_back(i);
//      visited.push_back(-1);
      };
   visited.push_back(-1);
   };

heap.push_back(nodelist[0]);
visited[nodelist[0]]=1;
nvisited++;

/* investigate connections by adding them to heap */
while(nvisited<(int)nodelist.size())
   {
   if(heap.size()==0)
     {
     valid=0;
     break;
     };

   id=heap.back();heap.pop_back();

   /* current node is id, investigate its connections, add to heap */

   for(i=0;i<(int)nodes.at(id).connections.size();i++)
      {
      conn_to=nodes.at(id).connections.at(i).to;

      if(segments[conn_to]==segid&&visited[conn_to]==-1) // if connected to other node within same segment and not visited, add to heap
        {
        heap.push_back(conn_to);
        nvisited++;
        visited[conn_to]=1;
        };
      };

   };


return valid;
}// end of CHANNEL::validateSegment()




/* grow and assign segments starting from a node that has a largest radii */
void FEATURE::segment_growSegment_start(ATOM_NETWORK *uc, int nodeid){

segments[nodeid]=nsegments; // assign new segment number to a nw segment
segment_max_r.push_back(nodes.at(nodeid).max_radius); //store maximum r in a segment
segment_growSegment_cont(uc,nodeid,nodeid);

nsegments++;
}// end of segment_growSegment()

/*
class CONN{
public:
  int from, to;       // IDs of nodes involved in connectioon
  double length;      // Length of edge
  double max_radius;  // Radius of largest spherical probe that can travel along edge
  DELTA_POS deltaPos; // Change in unit cell
*/

/* initnodeid is id of a nod that was a seed for this segment */
void FEATURE::segment_growSegment_cont(ATOM_NETWORK *uc,int initnodeid,int nodeid){
int i,j;
int conn_to;
int update_conn_flag; // segment connection updated or new
double r;
SEGCONN segcon;

/* this if makes sure that segment is grown only within the sphere marked by the shere located on the seed node */
if(uc->calcDistanceXYZ(nodes.at(initnodeid).x,nodes.at(initnodeid).y,nodes.at(initnodeid).z,nodes.at(nodeid).x,nodes.at(nodeid).y,nodes.at(nodeid).z)<SPHERE_SCALE*nodes.at(initnodeid).max_radius)
{

segments[nodeid]=nsegments;

r=nodes.at(nodeid).max_radius;

for(i=0;i<(int)nodes.at(nodeid).connections.size();i++)
   {
   conn_to=nodes.at(nodeid).connections.at(i).to;

   if(segments[conn_to]==-1) // investigate connection to nodes not yet assigned to any segments
      {

//      if(nodes.at(conn_to).max_radius<r&&(nodes.at(conn_to).max_radius/nodes.at(nodeid).connections.at(i).max_radius)<NODE_R_EDGE_R_RATIO)
         {
         segment_growSegment_cont(uc,initnodeid,conn_to);
         }
         ;
// 
/*
      else if(nodes.at(conn_to).max_radius==r&&
             ((nodes.at(conn_to).max_radius/nodes.at(nodeid).connections.at(i).max_radius)<1.05||
             (r/nodes.at(nodeid).connections.at(i).max_radius)<1.05)
             )
         {
         segment_growSegment_cont(conn_to);
         };
*/


      } // end if(segments[conn_to]==-1)
      else if(segments[conn_to]!=-1&&segments[conn_to]!=nsegments) // do if encounter segment-segment boundary
      {
      segcon.from=nodeid; // assign value to new segment connection
      segcon.to=conn_to;
      segcon.from_seg=nsegments;
      segcon.to_seg=segments[conn_to];
      segcon.max_radius=nodes.at(nodeid).connections.at(i).max_radius;
      segcon.length=nodes.at(nodeid).connections.at(i).length;
      segcon.merged=0;
      
      segment_conn_fulllist.push_back(segcon); // store every intersegment connection
                                               // this will be used to identify features
                                               // connecting different segments

      if(segment_conn.size()==0) {segment_conn.push_back(segcon);}
        else
           {
           update_conn_flag=0;
           for(j=0;j<(int)segment_conn.size();j++)
              {
//              if((segment_conn.at(j).from==nodeid&&segment_conn.at(j).to==conn_to)||(segment_conn.at(j).from==conn_to&&segment_conn.at(j).to==nodeid))
                if((segment_conn.at(j).from_seg==nsegments&&segment_conn.at(j).to_seg==segments[conn_to])
                      ||(segment_conn.at(j).from_seg==segments[conn_to]&&segment_conn.at(j).to==nsegments))
                 {
                 update_conn_flag=1;
                 if(segment_conn.at(j).max_radius<nodes.at(nodeid).connections.at(i).max_radius) 
                            segment_conn.at(j).max_radius=nodes.at(nodeid).connections.at(i).max_radius;
                 };
         /*        else
                 {             
                 segment_conn.push_back(segcon);
                 }; */
              };
           if(update_conn_flag==0) segment_conn.push_back(segcon);
           }; // end else if(segment_conn.size()==0)


      }; // end of else if
   };


}; // closes if(calculatedistance)....



}



/* segment detection algorithm identified only one (major) connection to other segment (with an edge with the largest corresponding r)
 * we need to post-process intersegment connecting nodes (ones with edges connecting to other segments) to obtain this info */

/* INFO
  vector <SEGCONN> segment_conn_fulllist; // lists all connections (edges) briging different segments
                                          // this will be used to determine topology of segments
  vector < vector<int> > segment_faces_conn;   // list of 'faces' for each segment, faces are clusters of edges that connect to another segment
  vector < vector<double> > segment_faces_max_r; // max r for each segment's face 
*/
void FEATURE::segment_identify_connections(){
int i,j,k;
int nfaces;
vector <int> segfaces;
vector <double> segfaces_max_r;

vector <int> nodes_in_segment; // list IDs of nodes assigned to segment
vector <int> face_nodes_in_segment; // list IDs of nodes assigned to segment and located on segment's "face", that is
                                    // nodes with connections to other segments
vector <int>  connected_segment_for_face_nodes_in_segment; // for each node, if it is a face node, store segment id to which it is connected
vector <SEGCONN> edges_from_segment; // list of edges that connect the current segment to other segments
SEGCONN segc; // inter-segment connection
vector <int> edges_in_face; // for each edge connecting to other semgent assign to a face, this svector stores face id for each edge

/* loop over all segments and identify all connections to other segments */

for(i=0;i<nsegments;i++)
   {
   nodes_in_segment.clear();
   edges_from_segment.clear();
   for(j=0;j<(int)nodes.size();j++)
      {
      if(segments[j]==i) nodes_in_segment.push_back(j);
      };
   face_nodes_in_segment.clear();
   connected_segment_for_face_nodes_in_segment.clear();
   face_nodes_in_segment.resize(nodes_in_segment.size(),0);
   connected_segment_for_face_nodes_in_segment.resize(nodes_in_segment.size(),-1);
   for(j=0;j<(int)nodes_in_segment.size();j++)
      {
      for(k=0;k<(int)nodes.at(nodes_in_segment[j]).connections.size();k++)
         {
         if(segments[nodes.at(nodes_in_segment[j]).connections.at(k).to]!=i) 
             {
             face_nodes_in_segment[j]++; 
             connected_segment_for_face_nodes_in_segment[j]=segments[nodes.at(nodes_in_segment[j]).connections.at(k).to];
        
             segc.from=nodes_in_segment[j]; // this is ID id in channel)
             segc.to=nodes.at(nodes_in_segment[j]).connections.at(k).to; // node id in channel
             segc.from_seg=i; 
             segc.to_seg=segments[nodes.at(nodes_in_segment[j]).connections.at(k).to];
             segc.max_radius=nodes.at(nodes_in_segment[j]).connections.at(k).max_radius;
             segc.deltaPos=nodes.at(nodes_in_segment[j]).connections.at(k).deltaPos;
             edges_from_segment.push_back(segc);
             };
         };
      };

// now connected_segment_for_face_nodes_in_segment has a list of segments each node in the current segment is conneed to (-1 if not) 
// we need now to find "clusters" of these nodes that will correspond to different connections;
// edges_from_segment stores all edges going out from a segment to other segments

   int nedges_tobe_assigned=edges_from_segment.size();
   edges_in_face.clear();
   edges_in_face.resize(nedges_tobe_assigned,-1);

   nfaces=0;
   int neigh_seg; // id of neighboring segment
//   while(nedges_tobe_assigned>-1)
     {
     for(k=0;k<(int)edges_from_segment.size();k++)
        {
        if(edges_in_face[k]==-1) // find first edge not assigned to face
          {
          neigh_seg=edges_from_segment[k].to_seg;
          edges_in_face[k]=nfaces;
          nedges_tobe_assigned--;
          segment_identify_connections_recurr(&nodes_in_segment,&edges_in_face,&nedges_tobe_assigned,nfaces,neigh_seg,&edges_from_segment,k);
          nfaces++;
          };
//        nfaces++;
        };
     };

/*
   cout << "Segment " << i <<": " << nodes_in_segment.size() << " nodes, " << edges_in_face.size() << " edges, " << nfaces << " faces, nodes-to-go " << nedges_tobe_assigned << "\n";
   for(k=0;k<(int)edges_from_segment.size();k++)
      {
      cout << "Edge " << edges_from_segment[k].from << "  " << edges_from_segment[k].to << "  Segments:" << edges_from_segment[k].from_seg << "  " << edges_from_segment[k].to_seg << "    "; edges_from_segment[k].deltaPos.print();
      };
*/


   }; // end of loop over all segments

} // end of segment_identify_connections()


/* segment_identify_connections_recurr() identifed edges that belong to particular face
 * */
void FEATURE::segment_identify_connections_recurr(vector <int> *nodeid,vector <int> *edgeinface,int *nodes_togo,int curr_face,int curr_neigh_seg, vector <SEGCONN> *edges,int curr_edge)
{
int i;
// check node at current edge, if not assigned to face, assign, consider its neighbors


for(i=0;i<(int)edges->size();i++)
   {
   if(i!=curr_edge&&edgeinface->at(i)==-1) // only consider if not yet assigned to face
     {
     //only consider pairs of edges that connect to the same segment 
     if(curr_neigh_seg==edges->at(i).to_seg&&edges->at(curr_edge).deltaPos.equals(edges->at(i).deltaPos))
       {
       // for two edges connecting to another segment, two edges are on the same face
       // only if starting nodes of this face overlap
//       int node1id=nodeid->at(edges->at(curr_edge).from);
//       int node2id=nodeid->at(edges->at(i).from);

//       int node1id=(edges->at(curr_edge).from);
//       int node2id=(edges->at(i).from);

//       if(nodes[node1id].max_radius+nodes[node2id].max_radius>calcEuclideanDistance(nodes[node1id].x,nodes[node1id].y,nodes[node1id].z,nodes[node2id].x,nodes[node2id].y,nodes[node2id].z))
          { // starting nodes of edges overlap, assign edge to the current face
          edgeinface->at(i)=curr_face;
          *nodes_togo=(*nodes_togo)--;
          segment_identify_connections_recurr(nodeid,edgeinface,nodes_togo,curr_face,curr_neigh_seg,edges,i);
          };


       }; // end of if - overlapping nodes

     }; // end of if - the same segments
   }; // end of loop over all edges



} // end of segment_identify_connections_recurr()




/* updated a vector with stores segment id for each node in voronoi network */
/* if flag==1 that it does the same but for features not segments */
void FEATURE::update_node_segmentinfo(vector <int> *segid,int *curr_id,int flag)
{
	int i;
	//loop over all nodes in channel system and update voronoi network nodes

	for(i=0;i<(int)nodes.size();i++)
	{
		if(flag==0)
		{
			segid->at(reverseIDMappings[i])=segments[i]+(*curr_id);
		}else{
			segid->at(reverseIDMappings[i])=feature_node_mapping[i]+(*curr_id);
		};
	};

	if(flag==0) {
		*curr_id+=nsegments;
	}else{
		*curr_id+=nfeatures;
	};

} // end of update_node_segmentinfo()



/* This function does perform "fake" merging - all segments are being copied as features
 * this is done to use all datastrucutres created for features to analyze segments  
*/
void FEATURE::merge_nothing()
{
	int i;

	nfeatures=nsegments;
	for(i=0;i<nsegments;i++)
	{
		features.push_back(i);
		feature_max_r.push_back(segment_max_r[i]);
	};

	for(i=0;i<(int)nodes.size();i++)
	{
		feature_node_mapping.push_back(segments[i]);
	};

	cout << "After merging nothing, nfeatures = " << nfeatures << "\n";

}



/* This is new merging that just mergest single node semgents */
void FEATURE::merge_singlenode_segments()
{
	int i,j,k;
	int mink=0;
	int count;

	int nsingle=0;
	double min=0;
	vector <int> nodelist;  // list of nodes in current segment
	vector <int> tobemerged; // flag for nodes to be merged with other segments
	vector <int> tobeextended; // flag for nodes that will be extended by toher nodes
	vector <int> listofsinglenodesegments;
	vector <int> tobemerged_segments;

	tobemerged.resize(nodes.size(),-1);
	tobeextended.resize(nodes.size(),0);
	tobemerged_segments.resize(nsegments,0);

	for(i=0;i<nsegments;i++)
	{
		nodelist.clear();
		for(j=0;j<(int)nodes.size();j++)
		{
			if(segments[j]==i) nodelist.push_back(j);
		};

		if(nodelist.size()==1)
		{


			if(tobeextended[nodelist[0]]==0)  // if the current node is not planned to be extended with other node, merge with others
			{
				tobemerged_segments[i]=1; // mark current segment as one to e merged

				for(k=0;k<(int)nodes.at(nodelist[0]).connections.size();k++) // iterate over all connections to find closes node
				{
					if(k==0) {min=nodes.at(nodelist[0]).connections.at(k).length;mink=nodes.at(nodelist[0]).connections.at(k).to;};
					if(k>0&&nodes.at(nodelist[0]).connections.at(k).length<min)
					{
						mink=nodes.at(nodelist[0]).connections.at(k).to;
						min=nodes.at(nodelist[0]).connections.at(k).length;
					};
				};
	// shortest connection leads to mink

				if(tobemerged[mink]>-1) mink=tobemerged[mink]; // if node to which we are going to merge the current node
													// also will be merged, merge into the other node

				tobemerged[nodelist[0]]=mink; // curr node to be connected to mink
				tobeextended[mink]=1;
				tobeextended[nodelist[0]]=1; // "write protect" current node

				nsingle++;
			};
		};

	};

/* now merging */

	nfeatures=nsegments-nsingle;

	feature_node_mapping.resize(nodes.size(),-1);

	count=0;
	for(i=0;i<nsegments;i++)
	{
		if(tobemerged_segments[i]==0)
		{
			features.push_back(count);
			feature_max_r.push_back(segment_max_r[i]);

			for(j=0;j<(int)nodes.size();j++)
			{
				if(segments[j]==i) feature_node_mapping[j]=count;
			};

			count++;
		};
	};
	count=0;
	for(i=0;i<(int)nodes.size();i++)
	{
		if(tobemerged[i]>-1)
		{
			feature_node_mapping[i]=feature_node_mapping[tobemerged[i]];
			count++;
		};
	};


	for(i=0;i<(int)nodes.size();i++)
	{
		if(feature_node_mapping[i]==-1) 
		{
			cout << "Feature construction failed. Abort."<< "\n";abort();
		}
	}

	cout << "After merging single segments, nfeatures = " << nfeatures << "\n";


} // end of merge_singlenode_segments()




int FEATURE::number_nodes_in_segment(int segid)
{
	int i;
	int count=0;

	for(i=0;i<(int)nodes.size();i++)
	{
		if(segments[i]==segid) count++;
	};

	return count;
}






/* merge_segments function goes over all identified segments and tries to merge some into larger blocks */

void FEATURE::merge_segments()
{
	int i;
	int current_segment;
	SEGCONN featcon;
	segment_merge_list.resize(nsegments,-1);

	nfeatures=0;


	/* loop over all segment connections and mark ones to be merge
	* Merging criteria - restriction is smaller than % of average max radii for each segment
	* */
	for(i=0;i<(int)segment_conn.size();i++)
	{
		if(segment_conn.at(i).max_radius>FEAT_NODE_R_EDGE_R_RATIO*0.5*(segment_max_r[(segment_conn.at(i).from_seg)]+segment_max_r[(segment_conn.at(i).to_seg)]))
			segment_conn.at(i).merged=1;
	/* new function to merge short edges: */

	//   if(segment_conn.at(i).length<max(segment_max_r[(segment_conn.at(i).from_seg)],segment_max_r[(segment_conn.at(i).to_seg)]))
	//       segment_conn.at(i).merged=1;

	};

	/* after connections of to-be-merged segments have been identified, we build
	* features, this includes 'flood fill' of segments */ 

	do{
		//find first segment not assigned to feature 
		current_segment=merge_findNotAssignedFragment();  
		//creature new feature
		features.push_back(nfeatures); // number of segments after merging 
		//grow
		//  cout << current_segment << "\n";
		merge_newSegmentGrow(current_segment);
		//increase counter (number of segments)
		nfeatures++;
	}while(merge_findNotAssignedFragment()!=-1);

	// update feature radii vector
	feature_max_r.resize(nfeatures,0);
	for(i=0;i<(int)segment_merge_list.size();i++)
	{
		if(feature_max_r[segment_merge_list[i]]<segment_max_r[i]) feature_max_r[segment_merge_list[i]]=segment_max_r[i]; 
	};
	// update feature information on all nodes in a channel
	for(i=0;i<(int)nodes.size();i++)
	{
		feature_node_mapping.push_back(segment_merge_list[segments[i]]);
	};

	/* having identified features we need to convert segments connectivity info
	* into feature connectivity */

	for(i=0;i<(int)segment_conn.size();i++)
	{
		if(segment_conn.at(i).merged!=1)
		{
			featcon.from_seg=segment_merge_list[segment_conn.at(i).from_seg];
			featcon.to_seg=segment_merge_list[segment_conn.at(i).to_seg];
			featcon.max_radius=segment_conn.at(i).max_radius;
			feature_conn.push_back(featcon);
		};
	};

} // end of CHANNEL:merge_segments()



void FEATURE::merge_newSegmentGrow(int segid)
{
	int k;
	int conn_to;

	segment_merge_list[segid]=nfeatures;


	for(k=0;k<(int)segment_conn.size();k++)
	{
		if((segment_conn.at(k).from_seg==segid||segment_conn.at(k).to_seg==segid)&&segment_conn.at(k).merged==1)
		{
			if(segment_conn.at(k).from_seg==segid) 
				{conn_to=segment_conn.at(k).to_seg;}
			else {conn_to=segment_conn.at(k).from_seg;};

			if(segment_merge_list.at(conn_to)==-1) merge_newSegmentGrow(conn_to);


		}; // end if
	}; //end for loop


}// end of merge_newSegmentGrow



int FEATURE::merge_findNotAssignedFragment()
{
	int i;
	int r=-1;
	for(i=segment_merge_list.size()-1;i>-1;i--)
	{
		if(segment_merge_list.at(i)==-1) r=i;
	};
	return r;
} //end merge_findNotAssignedFragment()


void FEATURE::segment_saveVis(fstream &output)
{
	int k;

	if(!output.is_open()){
		cerr << "Error: File stream needed to print segment information was not open." << "\n"
			<< "Exiting ..." << "\n";
		exit(1);
	}
	else{

		for(k=0;k<nsegments;k++)
		{
			output << "\n" << "Segment " << k << " with max_r of " << segment_max_r.at(k) << "\n";
	// Draw the components located in each unit cell
			for(unsigned int i = 0; i < unitCells.size(); i++){
				vector<int> nodeIDs = ucNodes.at(i);
				DELTA_POS disp = unitCells.at(i);

	// Iterate over all nodes in the unit cell
				for(unsigned int j = 0; j < nodeIDs.size(); j++){
					DIJKSTRA_NODE curNode = nodes.at(nodeIDs.at(j));

					if(segments.at(nodeIDs.at(j))==k)
					{
	// Find node coordinates
						double xCoord = curNode.x + v_a.x*disp.x + v_b.x*disp.y + v_c.x*disp.z;
						double yCoord = curNode.y + v_a.y*disp.x + v_b.y*disp.y + v_c.y*disp.z;
						double zCoord = curNode.z + v_a.z*disp.x + v_b.z*disp.y + v_c.z*disp.z;

	// Command used to draw node
//        output << "draw sphere {" << xCoord <<  " " << yCoord << " " << zCoord << "} radius " <<  curNode.max_radius << "\n";	
						output << "draw sphere {" << xCoord <<  " " << yCoord << " " << zCoord << "} radius " <<  0.1 << "\n";

					};//printing segment block


				} // ends for loop

			}
		}
	}
}


/** Calcluate volume of channel/segment/feature
 *  It uses reconstructed feature, puts it in a box
 *  and performs MC integration */
 
double FEATURE::calculateVolume(){
 
  long int i,j,counter=0;  
  double volume;

  XYZ max,min;
  XYZ box;
  XYZ pt;
  srand(randSeed2);
  SPHERE sph;
  vector <SPHERE> object; // a vector of spheres that define object (segemnt/feature/channel)

    // Draw the components located in each unit cell
    for(unsigned int i = 0; i < unitCells.size(); i++){
      vector<int> nodeIDs = ucNodes.at(i);
      DELTA_POS disp = unitCells.at(i);
      
      // Iterate over all nodes in the unit cell
      for(unsigned int j = 0; j < nodeIDs.size(); j++){
	DIJKSTRA_NODE curNode = nodes.at(nodeIDs.at(j));
	
	// Find node coordinates

        if(dimensionality<1) 
           { // if feature is not a channel, use reconstrcted feature
	   sph.x = curNode.x + v_a.x*disp.x + v_b.x*disp.y + v_c.x*disp.z;
	   sph.y = curNode.y + v_a.y*disp.x + v_b.y*disp.y + v_c.y*disp.z;
	   sph.z = curNode.z + v_a.z*disp.x + v_b.z*disp.y + v_c.z*disp.z;
           } else {
           sph.x = nodes.at(nodeIDs.at(j)).x;sph.y = nodes.at(nodeIDs.at(j)).y;sph.z = nodes.at(nodeIDs.at(j)).z;
//           cout << "draw sphere { " << sph.x << "  " << sph.y << "  " << sph.z << " } radius 0.1 " << "\n"; 
           };

	sph.r = SPHERE_SCALE*curNode.max_radius;

        object.push_back(sph);

        if(object.size()==1)
          { // first sphere added
          min.x=sph.x-sph.r;max.x=sph.x+sph.r;
          min.y=sph.y-sph.r;max.y=sph.y+sph.r;
          min.z=sph.z-sph.r;max.z=sph.z+sph.r;
          }else
          {
          if(min.x>(sph.x-sph.r)) min.x=sph.x-sph.r;
          if(max.x<(sph.x+sph.r)) max.x=sph.x+sph.r;
          if(min.y>(sph.y-sph.r)) min.y=sph.y-sph.r;
          if(max.y<(sph.y+sph.r)) max.y=sph.y+sph.r;
          if(min.z>(sph.z-sph.r)) min.z=sph.z-sph.r;
          if(max.z<(sph.z+sph.r)) max.z=sph.z+sph.r;
          };
		
      }
    } // end of analysis of object




  box.x=max.x-min.x;
  box.y=max.y-min.y;
  box.z=max.z-min.z;
 
  volume=box.x*box.y*box.z;

  cout << "Segment Box volume= " << volume << "\n"; 

// Now MC integration of volume

// example of random generation
// //     double theta = (rand()*1.0/RAND_MAX)*2*pi;

  for(i=0;i<N_MC_SAMPLE;i++)
     {

     pt.x = (rand()*1.0/RAND_MAX)*box.x+min.x;
     pt.y = (rand()*1.0/RAND_MAX)*box.y+min.y;
     pt.z = (rand()*1.0/RAND_MAX)*box.z+min.z;


     for(j=0;j<(int)object.size();j++)
        {
        if(calcEuclideanDistance(pt.x,pt.y,pt.z,object[j].x,object[j].y,object[j].z)<object[j].r)
           {
           counter++;
           break;
           };
        };

     };

volume=volume*((double)(N_MC_SAMPLE-counter))/(double)N_MC_SAMPLE;

return volume;

}
/* end of calculateVolume */

// new version from Thomas
// VORONOI_NETWORK *origVNet, DIJKSTRA_NETWORK *origDNet

void FEATURE::calculateVolumeForAllFeatures(VORONOI_NETWORK *origVNet,DIJKSTRA_NETWORK *origDNet){

    vector< vector<int> > nodesForFeatures = vector < vector<int> > (nfeatures, vector<int> ());
    vector<FEATURE> featChannels = vector<FEATURE> ();

    for(unsigned int i = 0; i < feature_node_mapping.size(); i++){
      nodesForFeatures[feature_node_mapping[i]].push_back(reverseIDMappings.find(i)->second);
    }

    int zeroVecs [3][3] = {{0, 0, 0},{ 0, 0, 0}, {0, 0, 0}};

    for(int i = 0; i < nfeatures; i++){
       DIJKSTRA_NETWORK featNet;
       DIJKSTRA_NETWORK::filterDnetEdges(nodesForFeatures[i], origVNet, &featNet);
       vector<bool> infoStorage; vector<CHANNEL> channels;
       findChannels(&featNet, &infoStorage, &channels);
       int dim;
       if(channels.size() != 0){
         dim = channels[0].dimensionality;
         featChannels.push_back(FEATURE(nodesForFeatures[i], origDNet, dim, channels[0].basis));
       }
       else{
         dim = 0;
         featChannels.push_back(FEATURE(nodesForFeatures[i], origDNet, 0, zeroVecs));
       }
//       featChannels[i].writeToVMD("feature", initIndex + i, output);
       feature_volume.push_back(featChannels[i].calculateVolume());
//       cout << "Feature volume = " << featChannels[i].calculateVolume() << "\n";
    };

  }






/** Calcluate surface area of channel/segment/feature
 *  It uses reconstructed feature, puts it in a box
 *  and performs MC integration */
 
double FEATURE::calculateSurfaceArea_OLD(ATOM_NETWORK *uc,vector<DIJKSTRA_NODE> *ChanNodeList,vector <int> *featnodemap,int curr_f,char *name,int fid){
 
  long int i,j,k,counter=0;  
  double volume;
  int step;
  double r_tr; // threshold for r
  XYZ max,min;
  XYZ box;
  XYZ pt;

  vector <double> SA_holo;  // stores hologram for each feature
  vector <double> SA_holo_otherf; // as above but excluding surface that is covered by other segments

  srand(randSeed2);
  SPHERE sph,sph2;
  vector <SPHERE> object; // a vector of spheres that define object (segemnt/feature/channel)

  vector <SPHERE> other_object; // vectors of spheres corresponding to other segments

    // Draw the components located in each unit cell
    for(unsigned int i = 0; i < unitCells.size(); i++){
      vector<int> nodeIDs = ucNodes.at(i);
      DELTA_POS disp = unitCells.at(i);
      
      // Iterate over all nodes in the unit cell
      for(unsigned int j = 0; j < nodeIDs.size(); j++){
	DIJKSTRA_NODE curNode = nodes.at(nodeIDs.at(j));
	
	// Find node coordinates

        if(dimensionality<1) 
           { // if feature is not a channel, use reconstrcted feature
	   sph.x = curNode.x + v_a.x*disp.x + v_b.x*disp.y + v_c.x*disp.z;
	   sph.y = curNode.y + v_a.y*disp.x + v_b.y*disp.y + v_c.y*disp.z;
	   sph.z = curNode.z + v_a.z*disp.x + v_b.z*disp.y + v_c.z*disp.z;
           } else {
           sph.x = nodes.at(nodeIDs.at(j)).x;sph.y = nodes.at(nodeIDs.at(j)).y;sph.z = nodes.at(nodeIDs.at(j)).z;
//           cout << "draw sphere { " << sph.x << "  " << sph.y << "  " << sph.z << " } radius 0.1 " << "\n"; 
           };

	sph.r = SPHERE_SCALE*curNode.max_radius;

        object.push_back(sph);

        if(object.size()==1)
          { // first sphere added
          min.x=sph.x-sph.r;max.x=sph.x+sph.r;
          min.y=sph.y-sph.r;max.y=sph.y+sph.r;
          min.z=sph.z-sph.r;max.z=sph.z+sph.r;
          }else
          {
          if(min.x>(sph.x-sph.r)) min.x=sph.x-sph.r;
          if(max.x<(sph.x+sph.r)) max.x=sph.x+sph.r;
          if(min.y>(sph.y-sph.r)) min.y=sph.y-sph.r;
          if(max.y<(sph.y+sph.r)) max.y=sph.y+sph.r;
          if(min.z>(sph.z-sph.r)) min.z=sph.z-sph.r;
          if(max.z<(sph.z+sph.r)) max.z=sph.z+sph.r;
          };
		
      }
    } // end of analysis of object



// Create a list of spheres that overlap with the current feature
/*
 for(j=0;j<(int)ChanNodeList->size();j++)
    {
    int count=0;
    int count2=0;
    Point node;
    node.x=ChanNodeList->at(j).x; node.y=ChanNodeList->at(j).y; node.z=ChanNodeList->at(j).z;
    node=uc->shiftXYZInUC(node);

    for(int i = 0; i < (int)object.size(); i++){
       Point p;
       p.x=object[i].x;p.y=object[i].y;p.z=object[i].z;
       Point p_shifted=uc->shiftXYZInUC(p);

       double dist = uc->calcDistanceXYZ(p.x,p.y,p.z,node.x,node.y,node.z);
       if(dist<0.01&&fabs(ChanNodeList->at(j).max_radius*SPHERE_SCALE-object[i].r)<0.01) 
          {
          count++; break;
          };

       if(dist<ChanNodeList->at(j).max_radius*SPHERE_SCALE+object[i].r)
         {
         count2++;
         };

       };

    if(count==0) // if node is not a part of the current feature
      {
      if(count2>0) {
        //origDNet->nodes[j] is an sphere overlapping with the current feature
        sph2.r=ChanNodeList->at(j).max_radius*SPHERE_SCALE;
        sph2.x=node.x; sph2.y=node.y; sph2.z=node.z;
        other_object.push_back(sph2);
        };
      };

    };
*/

// New version of code to find overlapping nodes

vector <int> conn_feat;

 for(j=0;j<(int)ChanNodeList->size();j++)
    {
    if(featnodemap->at(j)==curr_f)
      {
      for(i=0;i<(int)(ChanNodeList->at(j).connections.size());i++)      
         {
         if(featnodemap->at(ChanNodeList->at(j).connections.at(i).to)!=curr_f)
            {
            int flag=0;
            for(k=0;k<(int)conn_feat.size();k++)
               {
               if(conn_feat[k]==featnodemap->at(ChanNodeList->at(j).connections.at(i).to)) flag=1;
               };
            if(flag==0) conn_feat.push_back(featnodemap->at(ChanNodeList->at(j).connections.at(i).to));
            };
         };
      };
    }

for(i=0;i<(int)conn_feat.size();i++)
 for(j=0;j<(int)ChanNodeList->size();j++)
    {
    if(featnodemap->at(j)==i)
       {
       Point node;
       node[0]=ChanNodeList->at(j).x; node[1]=ChanNodeList->at(j).y; node[2]=ChanNodeList->at(j).z;
       node=uc->shiftXYZInUC(node);

       sph2.r=ChanNodeList->at(j).max_radius*SPHERE_SCALE;
       sph2.x=node[0]; sph2.y=node[1]; sph2.z=node[2];
       other_object.push_back(sph2);

       };
    };



cout << "Outside spheres list contain " << other_object.size() << " spheres (of " << ChanNodeList->size() << " nodes)." << "\n";


//

  box.x=max.x-min.x;
  box.y=max.y-min.y;
  box.z=max.z-min.z;
 
  volume=box.x*box.y*box.z;

//  cout << "Segment Box volume= " << volume << "\n"; 

// Now MC integration of surface

//int pos_origDNet;

SA_holo.resize(NBINS_SA_HOLO,0);
SA_holo_otherf.resize(NBINS_SA_HOLO,0);

for(step=0;step<NBINS_SA_HOLO;step++) // loop over all bins for surface area holograms
   {
   r_tr=(double)step*SA_HOLO_STEP;

   counter=0;

   for(j=0;j<(int)object.size();j++)
      {
      if(object[j].r>r_tr) { counter++;break;}; // at least one sphere larger than r_tr, sampling makes sense
      };


   if(counter>0)
     { // sampling makes sense

//     double totalSA = 0;
//     int sampleCount = 0;

     // Sample around each of the object's spheres
     for(int i = 0; i < (int)object.size(); i++){


        // now go to sampling

//        cout << "Object " << i << "\n"; cout.flush();

        int count = 0;
        int count_otherf=0;

        if(object[i].r-r_tr>0)
        for(int j = 0; j < N_MC_SA_SAMPLE; j++){

        bool overlap=false;
        bool overlap_otherf=false; // overlap with other features

        // Randomly sample point on a sphere of radius 1
        double theta = (rand()*1.0/RAND_MAX)*2*PI;
        double cosphi = 1.0 - (rand()*1.0/RAND_MAX)*2.0;
        double phi = acos(cosphi);

        // Convert spherical coordinates to xyz coordinates
        double xpoint = sin(phi)*cos(theta); //Don't need abs(sin(phi)) becase phi is from 0 to pi
        double ypoint = sin(phi)*sin(theta);
        double zpoint = cosphi;	
      
        // Scale the point to lie on a radius of r_probe + r_atom#i
        xpoint *= (object[i].r-r_tr);
        ypoint *= (object[i].r-r_tr);
        zpoint *= (object[i].r-r_tr);

        // move the point to the segment

        pt.x=object[i].x+xpoint;
        pt.y=object[i].y+ypoint;
        pt.z=object[i].z+zpoint;


        if(dimensionality<1)
          { // if the feature is not a channel, just check overlap with other spheres

          for(k=0;k<(int)object.size();k++)
             {
             if(k!=i)
                {
                if(calcEuclideanDistance(pt.x,pt.y,pt.z,object[k].x,object[k].y,object[k].z)<(object[k].r-r_tr)- 0.00001)
                   {
                   overlap=true; break;
                   };
                };
             };

          if(overlap==false) count++;


          //otherlap with other segments

         if(overlap==false)
         for(int jj=0;jj<(int)other_object.size();jj++)
            {
            Point p;
            p[0]=pt.x;p[1]=pt.y;p[2]=pt.z;
            Point p_shifted=uc->shiftXYZInUC(p);

            double dist = uc->calcDistanceXYZ(p[0],p[1],p[2],other_object[jj].x,other_object[jj].y,other_object[jj].z);

            if(dist<other_object[jj].r) {
                 overlap_otherf=true; break;
                 };

           };

 
          if(overlap_otherf==false) count_otherf++;

          }else{
          // sample points in case of channels, have to use pbc here

          // Convert (x,y,z) coordinates to ones relative to the unit cell vectors
          Point abc_coords = uc->xyz_to_abc(pt.x, pt.y, pt.z);

          abc_coords=uc->shiftABCInUC(abc_coords);
/*
          if(abc_coords.x<0||abc_coords.x>=1||abc_coords.y<0||abc_coords.y>=1||abc_coords.z<0||abc_coords.z>=1)
             {
             overlap=true;           
             };
*/
          if(overlap!=true)
            {
            for(k=0;k<(int)object.size();k++)
               {
               if(k!=i)
                  {
                  Point new_abc=uc->xyz_to_abc(object[k].x, object[k].y, object[k].z);
                  double centerDist = uc->calcDistanceABC(new_abc[0],new_abc[1],new_abc[2],abc_coords[0],abc_coords[1],abc_coords[2]);
                  if(centerDist < (object[k].r-r_tr) - 0.00001){
                     overlap=true; break;
                     };
                  };
               };

            };
          if(overlap==false) count++;

         // check otherlap with other segments
  
         if(overlap==false)
         for(int jk=0;jk<(int)other_object.size();jk++)
            {
            Point p;
            p[0]=pt.x;p[1]=pt.y;p[2]=pt.z;
            Point p_shifted=uc->shiftXYZInUC(p);

            double dist = uc->calcDistanceXYZ(p[0],p[1],p[2],other_object[jk].x,other_object[jk].y,other_object[jk].z);

            if(dist<other_object[jk].r) {
                 overlap_otherf=true; break;
                 };

           };
          

          if(overlap_otherf == false) 
	    count_otherf++;

          }; // end if, for channels (dim>0)

     
       }; // ends for(int j = 0; j < N_MC_SA_SAMPLE; j++){

//     cout << " count= " << count << "    cont_otherf= " << count_otherf << "\n";

     SA_holo[step]+= ((double)count/(double)N_MC_SA_SAMPLE)* 4.0 * PI * pow(object[i].r-r_tr,2);
     SA_holo_otherf[step]+= ((double)count_otherf/(double)N_MC_SA_SAMPLE)* 4.0 * PI * pow(object[i].r-r_tr,2);
     }; // ends for(int i = 0; i < (int)object.size() 

     }else   // from if(counter>0)
     {
     // surface area is 0
     SA_holo[step]=0;
     }; // ends  if(counter>0)

   }; // ends for(step=0;step<NBINS_SA_HOLO;step++)


cout << "##SA_holo:  " << name << "  Feature=  " << fid << "   ";
for(i=0;i<NBINS_SA_HOLO;i++)
   {
   cout << SA_holo[i] << "   ";
   };

cout << " Other:  ";
for(i=0;i<NBINS_SA_HOLO;i++)
   {
   cout << SA_holo_otherf[i] << "   ";
   };

cout << "\n";

//volume=volume*((double)(N_MC_SAMPLE-counter))/(double)N_MC_SAMPLE;

return SA_holo[0];


}
/* end of OLDcalculateSurfaceArea */



/* New version of feature surface area calcuation */

double FEATURE::calculateSurfaceArea(ATOM_NETWORK *uc,vector<DIJKSTRA_NODE> *ChanNodeList,vector <int> *featnodemap,int curr_f,char *name,int fid){
 
  long int i,j,k,counter=0;  
  int step;
  double r_tr; // threshold for r
  Point pt;

  vector <double> SA_holo;  // stores hologram for each feature
  vector <double> SA_holo_otherf; // as above but excluding surface that is covered by other segments

  srand(randSeed2);
  SPHERE sph,sph2;
  vector <SPHERE> object; // a vector of spheres that define object (segemnt/feature/channel)

  vector <SPHERE> other_object; // vectors of spheres corresponding to other segments



// Create object (current feature) and a list of spheres that overlap with the current feature

vector <int> conn_feat;

 for(j=0;j<(int)ChanNodeList->size();j++)
    {
    if(featnodemap->at(j)==curr_f)
      {
      /* if current feature, build object */
      sph.x=ChanNodeList->at(j).x; sph.y=ChanNodeList->at(j).y; sph.z=ChanNodeList->at(j).z;
      sph.r=SPHERE_SCALE*ChanNodeList->at(j).max_radius;

      object.push_back(sph);

      /* check feature two which it is connected */
      for(i=0;i<(int)(ChanNodeList->at(j).connections.size());i++)      
         {
         if(featnodemap->at(ChanNodeList->at(j).connections.at(i).to)!=curr_f)
            {
            int flag=0;
            for(k=0;k<(int)conn_feat.size();k++)
               {
               if(conn_feat[k]==featnodemap->at(ChanNodeList->at(j).connections.at(i).to)) flag=1;
               };
            if(flag==0) conn_feat.push_back(featnodemap->at(ChanNodeList->at(j).connections.at(i).to));
            };
         };
      };
    }

for(i=0;i<(int)conn_feat.size();i++)
 for(j=0;j<(int)ChanNodeList->size();j++)
    {
    if(featnodemap->at(j)==i)
       {
       sph2.x=ChanNodeList->at(j).x; sph2.y=ChanNodeList->at(j).y; sph2.z=ChanNodeList->at(j).z;
       sph2.r=ChanNodeList->at(j).max_radius*SPHERE_SCALE;

       int flag=0;
       for(k=0;k<(int)object.size();k++)
           {
           Point pp;
           pp[0]=ChanNodeList->at(j).x; pp[1]=ChanNodeList->at(j).y; pp[2]=ChanNodeList->at(j).z;
           Point pt_s=uc->minimizePointDistance(pp,object[k].x, object[k].y, object[k].z); // pt_s shifted
           double d = calcEuclideanDistance(pt_s[0],pt_s[1],pt_s[2],object[k].x, object[k].y, object[k].z);
           if(d<object[k].r+sph2.r) 
              {
              flag=1; break;
              };
           };
       if(flag==0) other_object.push_back(sph2);

       };
    };



cout << "Outside spheres list contain " << other_object.size() << " spheres (of " << ChanNodeList->size() << " nodes)." << "\n";



// Now MC integration of surface


SA_holo.resize(NBINS_SA_HOLO,0);
SA_holo_otherf.resize(NBINS_SA_HOLO,0);

for(step=0;step<NBINS_SA_HOLO;step++) // loop over all bins for surface area holograms
   {
   r_tr=(double)step*SA_HOLO_STEP;

   counter=0;

   for(j=0;j<(int)object.size();j++)
      {
      if(object[j].r>r_tr) { counter++;break;}; // at least one sphere larger than r_tr, sampling makes sense
      };


   if(counter>0)
     { // sampling makes sense


     // Sample around each of the object's spheres
     for(int i = 0; i < (int)object.size(); i++){

        int count = 0;
        int count_otherf=0;

        if(object[i].r-r_tr>0)
        for(int j = 0; j < N_MC_SA_SAMPLE; j++){

        bool overlap=false;
        bool overlap_otherf=false; // overlap with other features

        // Randomly sample point on a sphere of radius 1
        double theta = (rand()*1.0/RAND_MAX)*2*PI;
        double cosphi = 1.0 - (rand()*1.0/RAND_MAX)*2.0;
        double phi = acos(cosphi);

        // Convert spherical coordinates to xyz coordinates
        double xpoint = sin(phi)*cos(theta); //Don't need abs(sin(phi)) becase phi is from 0 to pi
        double ypoint = sin(phi)*sin(theta);
        double zpoint = cosphi;	
      
        // Scale the point to lie on a radius of r_probe + r_atom#i
        xpoint *= (object[i].r-r_tr);
        ypoint *= (object[i].r-r_tr);
        zpoint *= (object[i].r-r_tr);

        // move the point to the segment

        pt[0]=object[i].x+xpoint;
        pt[1]=object[i].y+ypoint;
        pt[2]=object[i].z+zpoint;

//             if(fid==0) cout << "### draw point { "<< pt.x << "  " << pt.y << "  " << pt.z<< " } " << "\n";
 
        Point newpt=uc->xyz_to_abc(pt);
        newpt=uc->shiftABCInUC(newpt);
        pt=uc->abc_to_xyz(newpt);
//        pt = uc->shiftXYZInUC(pt);
//             if(fid==0) cout << "### draw point { "<< pt.x << "  " << pt.y << "  " << pt.z<< " } " << "\n";

        for(k=0;k<(int)object.size();k++)
           {
           if(k!=i)
             {
//             Point pt_s=minimizePointDistance(pt,object[k].x, object[k].y, object[k].z, uc); // pt_s shifted
//             double centerDist = calcEuclideanDistance(pt_s.x,pt_s.y,pt_s.z,object[k].x, object[k].y, object[k].z); 

             double centerDist = uc->calcDistanceXYZ(pt[0],pt[1],pt[2],object[k].x, object[k].y, object[k].z);

             if(centerDist < (object[k].r-r_tr) - 0.00001){
                overlap=true; break;
                };
             };
           };

          if(overlap==false) count++;

         // check otherlap with other segments
  
         if(overlap==false)
         for(int jk=0;jk<(int)other_object.size();jk++)
            {

//           Point pt_s=minimizePointDistance(pt,other_object[jk].x,other_object[jk].y,other_object[jk].z,uc);
//            double dist = calcEuclideanDistance(pt_s.x,pt_s.y,pt_s.z,other_object[jk].x,other_object[jk].y,other_object[jk].z);

            double dist = uc->calcDistanceXYZ(pt[0],pt[1],pt[2],object[jk].x, object[jk].y, object[jk].z);

            if(dist<other_object[jk].r) {
                 overlap_otherf=true; break;
                 };

           };
          

          if(overlap==false&&overlap_otherf==false) 
             {
             count_otherf++;

//             if(fid==0) cout << "### draw point { "<< pt.x << "  " << pt.y << "  " << pt.z<< " } " << "\n"; 

             };


     
       }; // ends for(int j = 0; j < N_MC_SA_SAMPLE; j++){

     cout << " count= " << count << "    cont_otherf= " << count_otherf << "\n";

     SA_holo[step]+= ((double)count/(double)N_MC_SA_SAMPLE)* 4.0 * PI * pow(object[i].r-r_tr,2);
     SA_holo_otherf[step]+= ((double)count_otherf/(double)N_MC_SA_SAMPLE)* 4.0 * PI * pow(object[i].r-r_tr,2);
     }; // ends for(int i = 0; i < (int)object.size() 

     }else   // from if(counter>0)
     {
     // surface area is 0
     SA_holo[step]=0;
     }; // ends  if(counter>0)

   }; // ends for(step=0;step<NBINS_SA_HOLO;step++)


cout << "##SA_holo:  " << name << "  Feature=  " << fid << "   ";
for(i=0;i<NBINS_SA_HOLO;i++)
   {
   cout << SA_holo[i] << "   ";
   };

cout << " Other:  ";
for(i=0;i<NBINS_SA_HOLO;i++)
   {
   cout << SA_holo_otherf[i] << "   ";
   };

cout << "\n";

//volume=volume*((double)(N_MC_SAMPLE-counter))/(double)N_MC_SAMPLE;

return SA_holo[0];


}
/* end of calculateSurfaceArea */
