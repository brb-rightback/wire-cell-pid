#include "WCPPID/ProtectOverClustering.h"

using namespace WCP;

void WCPPID::Protect_Over_Clustering(std::vector<int>& to_be_checked, WCPPID::PR3DClusterSelection& live_clusters, std::map<WCPPID::PR3DCluster*, int>& map_cluster_parent_id,  std::map<int, std::vector<WCPPID::PR3DCluster*> >& map_parentid_clusters, WCP::ToyCTPointCloud& ct_point_cloud ){
  // update live_clusters,  map_parentid_clusters, map_cluster_parent_id;
  // keep the original id for the main cluster ...

  std::set<int> used_cluster_ids;
  for (auto it = live_clusters.begin(); it!=live_clusters.end(); it++){
    used_cluster_ids.insert((*it)->get_cluster_id());
    //std::cout << "haha: " << (*it)->get_cluster_id() << std::endl;
  }
  
  WCPPID::PR3DClusterSelection temp_live_clusters;
  int acc_cluster_id = 0;
  // std::cout << get_next_cluster_id(acc_cluster_id, used_cluster_ids) << std::endl;
  
  std::set<WCPPID::PR3DCluster*> examined_clusters;
  for (size_t i=0;i!=to_be_checked.size(); i++){
    int curr_main_cluster_id = to_be_checked.at(i);
    auto it = map_parentid_clusters.find(curr_main_cluster_id);
    if (it == map_parentid_clusters.end()) continue;
    for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
      WCPPID::PR3DCluster *temp_cluster = *it1;
      examined_clusters.insert(temp_cluster);
      map_cluster_parent_id.erase(temp_cluster);
      used_cluster_ids.erase(temp_cluster->get_cluster_id());
    }
    used_cluster_ids.insert(curr_main_cluster_id);
    map_parentid_clusters[curr_main_cluster_id].clear();
    
    // examine these clusters ...
    for (auto it1 = examined_clusters.begin(); it1 != examined_clusters.end(); it1++){
      WCPPID::PR3DCluster *temp_cluster = *it1;
      std::vector<WCP::SMGCSelection> vec_mcells;
      // hack for now ...
      //      if (temp_cluster->get_cluster_id() == curr_main_cluster_id)
	vec_mcells = temp_cluster->Examine_graph(ct_point_cloud);

      // check mcell size ...
      // int nsum = 0;
      // for (size_t j=0;j!=vec_mcells.size();j++){
      // 	nsum += vec_mcells.at(j).size();
      // }
      // std::cout << nsum << " " << temp_cluster->get_num_mcells() << std::endl;
      
      if (temp_cluster->get_cluster_id()==curr_main_cluster_id){
	int max_number_cells = vec_mcells.front().size();
	int main_id = 0;
	for (size_t j=1;j<vec_mcells.size();j++){
	  if (vec_mcells.at(j).size() > max_number_cells){
	    main_id = j;
	    max_number_cells = vec_mcells.at(j).size();
	  }
	}
	// test...
	// for (size_t j=0;j<vec_mcells.size();j++){
	//   double sum_q = 0;
	//   for (size_t k=0;k!=vec_mcells.at(j).size();k++){
	//     sum_q += vec_mcells.at(j).at(k)->get_q();
	//   }
	//   std::cout << j << " " << vec_mcells.at(j).size() << " " << sum_q << std::endl;
	// }
	
	// main cluster replacement 
	WCPPID::PR3DCluster *new_cluster = new PR3DCluster(curr_main_cluster_id);
	for (auto it2 = vec_mcells.at(main_id).begin(); it2!=vec_mcells.at(main_id).end();it2++){
	  new_cluster->AddCell(*it2, (*it2)->GetTimeSlice());
	}
	temp_live_clusters.push_back(new_cluster);
	map_cluster_parent_id[new_cluster] = curr_main_cluster_id;
	map_parentid_clusters[curr_main_cluster_id].push_back(new_cluster);

	// rest ...
	for (size_t j = 0; j!=vec_mcells.size(); j++){
	  if (j==main_id) continue;
	  acc_cluster_id = get_next_cluster_id(acc_cluster_id, used_cluster_ids);
	  WCPPID::PR3DCluster *new_cluster = new PR3DCluster(acc_cluster_id);
	  for (auto it2 = vec_mcells.at(j).begin(); it2!=vec_mcells.at(j).end();it2++){
	    new_cluster->AddCell(*it2, (*it2)->GetTimeSlice());
	  }
	  used_cluster_ids.insert(acc_cluster_id);
	  temp_live_clusters.push_back(new_cluster);
	  map_cluster_parent_id[new_cluster] = curr_main_cluster_id;
	  map_parentid_clusters[curr_main_cluster_id].push_back(new_cluster);
	}
	
	
      }else{
	for (size_t j = 0; j!=vec_mcells.size(); j++){
	  acc_cluster_id = get_next_cluster_id(acc_cluster_id, used_cluster_ids);
	  WCPPID::PR3DCluster *new_cluster = new PR3DCluster(acc_cluster_id);
	  for (auto it2 = vec_mcells.at(j).begin(); it2!=vec_mcells.at(j).end();it2++){
	    new_cluster->AddCell(*it2, (*it2)->GetTimeSlice());
	  }
	  used_cluster_ids.insert(acc_cluster_id);
	  temp_live_clusters.push_back(new_cluster);
	  map_cluster_parent_id[new_cluster] = curr_main_cluster_id;
	  map_parentid_clusters[curr_main_cluster_id].push_back(new_cluster);
	}
      }
    }
  }
   
  for (auto it = live_clusters.begin(); it != live_clusters.end(); it++){
    if (examined_clusters.find(*it)==examined_clusters.end())
      temp_live_clusters.push_back(*it);
  }
  
  for (auto it = examined_clusters.begin(); it!=examined_clusters.end(); it++){
    delete (*it);
  }
  
  live_clusters = temp_live_clusters;
  
   
}



int WCPPID::get_next_cluster_id(int acc_cluster_id, std::set<int>& used_cluster_ids){
  acc_cluster_id ++;
  while (used_cluster_ids.find(acc_cluster_id)!=used_cluster_ids.end()) acc_cluster_id++;
  return acc_cluster_id;
}
