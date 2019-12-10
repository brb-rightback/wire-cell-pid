 
std::tuple<int, WCPPID::PR3DCluster*, WCP::Opflash*> WCPPID::ToyFiducial::cosmic_tagger(WCP::OpflashSelection& flashes, WCPPID::PR3DCluster* main_cluster, std::vector<WCPPID::PR3DCluster*> additional_clusters, WCP::Opflash* main_flash, std::tuple<int, double, double, int>& bundle_info, WCP::Photon_Library *pl, int time_offset, int nrebin, float unit_dis, WCP::ToyCTPointCloud& ct_point_cloud, int run_no, int subrun_no, int event_no, bool flag_data, bool debug_tagger){


//	std::cout << "starting cosmic tagger ================================================================================" << std::endl;

	// 0 for original match results
	// 1 for TGM ...
	// 2 for STM ...
	int flag_background = 0;

	// bundle information ...
	int bundle_type = std::get<0>(bundle_info);
	double bundle_ks = std::get<1>(bundle_info);
	double bundle_chi2 = std::get<2>(bundle_info);
	int bundle_ndf = std::get<3>(bundle_info);
	//

	//Take in flash and cluster info
	std::vector<std::vector<WCPointCloud<double>::WCPoint>> extreme_points = main_cluster->get_extreme_wcps();
	double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;
	double main_offset_x = (main_flash->get_time() - time_offset)*2./nrebin*time_slice_width;

	//set tolerances
	double bundle_ks_stm_tol =     0.08;					//This is the tolerance below which a flash is considered a very good match and should not be considered
	double bundle_ks_tgm_tol =     0.00;					//TGM-with-flash is pure already, no need to screen by KS
	double bundle_ks_tgm_nof_tol = 0.08;
	double stm_flash_tol = 3.0*units::m;
	double tgm_flash_tol = 3.0*units::m;
	double stm_pe_frac_tol = 2.0;
	double tgm_pe_frac_tol = 2.0;
	std::vector<double> stm_tol_vec =     {0.0, 1.0, 1.0, 1.5, 1.0};	//x_ano, x_cat, ybot, ytop, z
	std::vector<double> tgm_tol_vec =     {2.0, 3.0, 3.0, 3.0, 3.0};
	std::vector<double> tgm_nof_tol_vec = {1.0, 1.5, 1.5, 1.5, 1.5};
	for(int i=0;i<5;i++){
		stm_tol_vec[i] *= units::cm;
		tgm_tol_vec[i] *= units::cm;
		tgm_nof_tol_vec[i] *= units::cm;
	}

	//PMT map and geometry info
	std::vector<int> OpDet_to_OpChannel_map = {4,2,5,0,6,1,3,9,11,7,12,8, 10,15,13,17,18,14,16,21,23,19,20,24,22,28,25,30,31,26,27,29};
	std::vector<int> OpChannel_to_OpDet_map = {3,5,1,6,0,2,4,9,11,7,12,8,10,14,17,13,18,15,16,21,22,19,24,20,23,26,29,30,25,31,27,28};
	std::vector<std::vector<double>> opDet_xyz_map = {{2.645, -28.625, 990.356},{2.682,  27.607, 989.712},{2.324, -56.514, 951.865},{2.458,  55.313, 951.861},{2.041, -56.309, 911.939},{2.265,  55.822, 911.066},{1.923, -0.722,  865.598},{1.795, -0.502,  796.208},{1.495, -56.284, 751.905},{1.559,  55.625, 751.884},{1.487, -56.408, 711.274},{1.438,  55.800, 711.073},{1.475, -0.051,  664.203},{1.448, -0.549,  585.284},{1.226,  55.822, 540.929},{1.479, -56.205, 540.616},{1.505, -56.323, 500.221},{1.116,  55.771, 500.134},{1.481, -0.875,  453.096},{1.014, -0.706,  373.839},{1.451, -57.022, 328.341},{0.913,  54.693, 328.212},{0.682,  54.646, 287.976},{1.092, -56.261, 287.639},{0.949, -0.829,  242.014},{0.658, -0.303,  173.743},{0.703,  55.249, 128.355},{0.821, -56.203, 128.179},{0.862, -56.615, 87.8695},{0.558,  55.249, 87.7605},{0.665,  27.431, 51.1015},{0.947, -28.576, 50.4745}};

	//Vectors that store STM / TGM / TGM-noflash cases
	std::vector<int> n_boundary_list;
	std::vector<Opflash*> candidate_flash_list;
	std::vector<double> candidate_offset_x_list;

	//Check whether the track is at least 15 cm long
	Point p0(extreme_points[0][0].x,extreme_points[0][0].y,extreme_points[0][0].z);
	Point p1(extreme_points[1][0].x,extreme_points[1][0].y,extreme_points[1][0].z);
	double distance = sqrt(pow(p0.x-p1.x,2)+pow(p0.y-p1.y,2)+pow(p0.z-p1.z,2));
	if(distance > 10*units::cm){
		//iterate over all the flashes to check each flash for STM/TGM conditions
		for (auto it = flashes.begin(); it!= flashes.end(); it++){
			Opflash *flash = (*it);
			double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;

			//Skip the already matched flash
			if(flash->get_time() == main_flash->get_time()){
				continue;
			}
			//Skip cases where the flash is clearly out of the beam window
			if(!inside_x_region(extreme_points,offset_x,0-tgm_tol_vec[0],256+tgm_tol_vec[1])){
				continue;
			}

			//Flash PE Prediction
			bool flag_good_bundle;
//			std::vector<double>& pred_pmt_light = main_bundle->get_pred_pmt_light();
//			WireCell2dToy::calculate_pred_pe(run_no, time_offset, nrebin, time_slice_width, pl, new_bundle, &pred_pmt_light, &additional_clusters, &other_clusters, &more_clusters, flag_good_bundle, flag_data);
			std::vector<double> pred_pmt_light = calculate_pred_pe(run_no, offset_x, pl, main_cluster, additional_clusters, flash, flag_data);

			//Compute the PE centroids for the flash PE and predicted PE
			double pred_pe_tot = 0;
			double flash_pe_tot = 0;
			double pred_pe_z_centroid = 0;
			double flash_pe_z_centroid = 0;
			for(int i=0;i<int(pred_pmt_light.size());i++){
				double flash_pe = flash->get_PE(i);
				flash_pe_tot += flash_pe;
				pred_pe_tot += pred_pmt_light[i];
				flash_pe_z_centroid += flash_pe*opDet_xyz_map[i][2]*units::cm;
				pred_pe_z_centroid += pred_pmt_light[i]*opDet_xyz_map[i][2]*units::cm;
			}
			flash_pe_z_centroid /= flash_pe_tot;
			pred_pe_z_centroid  /= pred_pe_tot;

			//Geometry evaluation
			//Call the helper function check_boundary to get the number of boundary intersections for a given flash
			//Only store a STM / TGM if the flash/predicted-PE roughly match
			int boundary_num_tgm = 0;
			int boundary_num_stm = 0;
			int boundary_num = 0;
			//TGM-with-flash case
			if(bundle_ks > bundle_ks_tgm_tol && flash_pe_tot < pred_pe_tot*tgm_pe_frac_tol && flash_pe_z_centroid > pred_pe_z_centroid-tgm_flash_tol && flash_pe_z_centroid < pred_pe_z_centroid+tgm_flash_tol){
				boundary_num_tgm = check_boundary(extreme_points, offset_x, &tgm_tol_vec);
			}
			//STM-with-flash case
			if(inside_x_region(extreme_points,offset_x,128-stm_tol_vec[0],256+stm_tol_vec[1]) && bundle_ks > bundle_ks_stm_tol && flash_pe_tot < pred_pe_tot*stm_pe_frac_tol && flash_pe_z_centroid > pred_pe_z_centroid-stm_flash_tol && flash_pe_z_centroid < pred_pe_z_centroid+stm_flash_tol){
				boundary_num_stm = check_boundary(extreme_points, offset_x, &stm_tol_vec);
			}
			//TGM gets priority over STM
			if(boundary_num_tgm==2)					{boundary_num =  2;}
			else if(boundary_num_stm==1)				{boundary_num =  1;}
			else if(boundary_num_tgm<0 || boundary_num_stm<0)	{boundary_num = -1;}

			//Store a STM / TGM tag if appropriate
			bool tgm = check_tgm(main_cluster,flash,offset_x,ct_point_cloud);
			if(boundary_num==1 || boundary_num==2){
				n_boundary_list.push_back(boundary_num);
				candidate_flash_list.push_back(flash);
			}
		}

		//Check whether a TGM w/o flash hypothesis works
		//Skip if the existing neutrino match is very good
		if(bundle_ks > bundle_ks_tgm_nof_tol){
			double offset_x = (main_flash->get_time() - time_offset)*2./nrebin*time_slice_width;
			double tx     = tgm_nof_tol_vec[0];
			double ty_bot = tgm_nof_tol_vec[1];
			double ty_top = tgm_nof_tol_vec[2];
			double tz     = tgm_nof_tol_vec[3];
			double max_offset =       1000000*units::cm;		//The max offset allowed in the +x direction.  The starting value is meant to be immediately overriden.
			double max_intersection = 1000000*units::cm;
			double min_intersection = 1000000*units::cm;
			double max_x =           -1000000*units::cm;
			double min_x =            1000000*units::cm;
			double pca0_x_min =       1000000*units::cm;
			double pca1_x_min =       1000000*units::cm;
			double max_x_anode =      0;
			double min_x_anode =      0;
			int pca_x_min_index =    -1;
			int pca_x_max_index =    -1;
			int pca_x_min_subindex = -1;
			for(int i=0;i<int(extreme_points.size());i++){
				for(int ii=0;ii<int(extreme_points[i].size());ii++){
					Point p(extreme_points[i][ii].x,extreme_points[i][ii].y,extreme_points[i][ii].z);
					if(i==0)     {pca0_x_min = std::min(pca0_x_min,(double) p.x);}
					else if(i==1){pca1_x_min = std::min(pca1_x_min,(double) p.x);}
				}
			}
			if(pca0_x_min < pca1_x_min){pca_x_min_index = 0;pca_x_max_index = 1;}	//only set the pca_x_max with the opposing pca end, so 0->1 and 1->0
			else			   {pca_x_min_index = 1;pca_x_max_index = 0;}

			//Iterate extreme points.  For each, compute the drift distance until the boundary+tolerance is crossed.
			//Set the max_offset variable to the lowest drift distance that causes any extreme point to (nearly) leave the boundary+tolerance.
			for(int i=0;i<int(extreme_points.size());i++){
				for(int ii=0;ii<int(extreme_points[i].size());ii++){
					Point p(extreme_points[i][ii].x,extreme_points[i][ii].y,extreme_points[i][ii].z);

					std::vector<double> boundary_points_xy_x = get_boundary_SCB_xy_x(p);
					std::vector<double> boundary_points_xy_y = get_boundary_SCB_xy_y(p);
					std::vector<double> boundary_points_xz_x = get_boundary_SCB_xz_x(p);
					std::vector<double> boundary_points_xz_z = get_boundary_SCB_xz_z(p);
					//Don't let the offset be larger than the distance to the cathode side
					max_offset = std::min(max_offset,boundary_points_xy_x[2]-p.x+tx);
					max_offset = std::min(max_offset,boundary_points_xz_x[2]-p.x+tx);
					max_intersection = std::min(max_intersection,boundary_points_xy_x[2]+tx);
					max_intersection = std::min(max_intersection,boundary_points_xz_x[2]+tx);
					min_intersection = std::min(min_intersection,boundary_points_xy_x[2]-tx);
					min_intersection = std::min(min_intersection,boundary_points_xz_x[2]-tx);
					max_x = std::max(max_x,(double) p.x);
					min_x = std::min(min_x,(double) p.x);
					if(min_x== (double) p.x){pca_x_min_subindex = ii;}
					max_x_anode = boundary_points_xy_x[0]+tx-0.01*units::cm;
					min_x_anode = boundary_points_xy_x[0]-tx+0.01*units::cm;
					//Conditions in Y and Z for whether the endpoint will drift into the SCB before the cathode side
					if     (p.y > boundary_points_xy_y[1]-ty_bot && p.y < boundary_points_xy_y[2]-ty_bot){
						double x0 = boundary_points_xy_x[1];	double y0 = boundary_points_xy_y[1]-ty_bot;
						double x1 = boundary_points_xy_x[2]+tx;	double y1 = boundary_points_xy_y[2]-ty_bot;
						double x_intersection = x0 + (x1-x0)*(p.y-y0)/(y1-y0);
						max_offset = std::min(max_offset,x_intersection-p.x);
						if(i==pca_x_max_index){
							max_intersection = std::min(max_intersection,x_intersection);
							x0 = boundary_points_xy_x[1];		y0 = boundary_points_xy_y[1]+ty_bot;
							x1 = boundary_points_xy_x[2]-tx;	y1 = boundary_points_xy_y[2]+ty_bot;
							x_intersection = x0 + (x1-x0)*(p.y-y0)/(y1-y0);
							min_intersection = std::min(min_intersection,x_intersection);
						}
					} else if(p.y > boundary_points_xy_y[3]+ty_top && p.y < boundary_points_xy_y[4]+ty_top){
						double x0 = boundary_points_xy_x[3]+tx;	double y0 = boundary_points_xy_y[3]+ty_top;
						double x1 = boundary_points_xy_x[4];	double y1 = boundary_points_xy_y[4]+ty_top;
						double x_intersection = x0 + (x1-x0)*(p.y-y0)/(y1-y0);
						max_offset = std::min(max_offset,x_intersection-p.x);
						if(i==pca_x_max_index){
							max_intersection = std::min(max_intersection,x_intersection);
							x0 = boundary_points_xy_x[3]-tx;	y0 = boundary_points_xy_y[3]-ty_top;
							x1 = boundary_points_xy_x[4];		y1 = boundary_points_xy_y[4]-ty_top;
							x_intersection = x0 + (x1-x0)*(p.y-y0)/(y1-y0);
							min_intersection = std::min(min_intersection,x_intersection);
						}
					}
					if     (p.z > boundary_points_xz_z[1]-tz && p.z < boundary_points_xz_z[2]-tz){
						double x0 = boundary_points_xz_x[1];	double z0 = boundary_points_xz_z[1]-tz;
						double x1 = boundary_points_xz_x[2]+tx;	double z1 = boundary_points_xz_z[2]-tz;
						double x_intersection = x0 + (x1-x0)*(p.z-z0)/(z1-z0);
						max_offset = std::min(max_offset,x_intersection-p.x);
						if(i==pca_x_max_index){
							max_intersection = std::min(max_intersection,x_intersection);
							x0 = boundary_points_xz_x[1];		z0 = boundary_points_xz_z[1]+tz;
							x1 = boundary_points_xz_x[2]-tx;	z1 = boundary_points_xz_z[2]+tz;
							x_intersection = x0 + (x1-x0)*(p.z-z0)/(z1-z0);
							min_intersection = std::min(min_intersection,x_intersection);
						}
					} else if(p.z > boundary_points_xz_z[3]+tz && p.z < boundary_points_xz_z[4]+tz){
						double x0 = boundary_points_xz_x[3]+tx;	double z0 = boundary_points_xz_z[3]+tz;
						double x1 = boundary_points_xz_x[4];	double z1 = boundary_points_xz_z[4]+tz;
						double x_intersection = x0 + (x1-x0)*(p.z-z0)/(z1-z0);
						max_offset = std::min(max_offset,x_intersection-p.x);
						if(i==pca_x_max_index){
							max_intersection = std::min(max_intersection,x_intersection);
							x0 = boundary_points_xz_x[3]-tx;	z0 = boundary_points_xz_z[3]-tz;
							x1 = boundary_points_xz_x[4];		z1 = boundary_points_xz_z[4]-tz;
							x_intersection = x0 + (x1-x0)*(p.z-z0)/(z1-z0);
							min_intersection = std::min(min_intersection,x_intersection);
						}
					}
				}
			}

			//Keep the offset within the boundary, and then check for TGM conditions.  If so, insert to the list, and use n_boudnary id of 3 to clarify that no flash was used.
			max_offset -= 0.01*units::cm;
			max_intersection -= 0.01*units::cm;
			min_intersection += 0.01*units::cm;
			//check cathode side
			int nboundaries = check_boundary(extreme_points,offset_x-max_offset, &tgm_nof_tol_vec);
			//check anode side
			double x_length = max_x-min_x;
			double max_x_anode_intersection = max_x_anode+x_length;
			double min_x_anode_intersection = min_x_anode+x_length;
			max_x_anode_intersection = std::min(max_x_anode_intersection,max_intersection);
			min_x_anode_intersection = std::max(min_x_anode_intersection,min_intersection);
			if(max_x_anode_intersection > min_x_anode_intersection){
				double adjusted_anode_middle = ((max_x_anode_intersection+min_x_anode_intersection)/2-x_length);
				double anode_boundary_offset = extreme_points[pca_x_min_index][pca_x_min_subindex].x-((max_x_anode_intersection+min_x_anode_intersection)/2-x_length);
				nboundaries = std::max(nboundaries,check_boundary(extreme_points,anode_boundary_offset, &tgm_nof_tol_vec));
			}

			if(nboundaries==2){
				n_boundary_list.push_back(3);
				candidate_offset_x_list.push_back(-offset_x+max_offset);
			}
		}

	}

	//Iterate the list and pick the highest priority tagging (TGM > TGM-noflash > STM > no tag)
	int tag_type = 0;
	WCP::Opflash* new_flash = main_flash;
	for(int n_match = 0; n_match<int(n_boundary_list.size());n_match++){
		if(n_boundary_list[n_match]==2)							{tag_type = 2; new_flash = candidate_flash_list[n_match];}
		else if(n_boundary_list[n_match]==3 && tag_type != 2)				{tag_type = 3; new_flash = NULL;}
		else if(n_boundary_list[n_match]==1 && tag_type != 2 && tag_type != 3)		{tag_type = 1; new_flash = candidate_flash_list[n_match];}
	}

	//debug_tagger = true;
	//Write to file for debugging purposes
	if(debug_tagger){
		for(int n_match = 0; n_match<int(n_boundary_list.size());n_match++){
			if(n_match < int(candidate_flash_list.size())){
				std::cout << "match type = " << n_boundary_list[n_match] << ", flash # = " << candidate_flash_list[n_match]->get_flash_id() << std::endl;		
			} else {
				std::cout << "match type = " << n_boundary_list[n_match] << std::endl;
			}
		}
		srand(time(NULL));
		int random_number = std::rand();
		std::ofstream debug_file;
		std::string path = "/uboone/data/users/lcoopert/cosmic_tagger/data";
//		debug_file.open(path+"/temp_numu_cc_all/numu_cc_tagger_results_"+std::to_string(random_number)+".txt", std::ios_base::app);
//		debug_file.open(path+"/temp_extbnb_all/extbnb_tagger_results_"+std::to_string(random_number)+".txt", std::ios_base::app);

//		debug_file.open(path+"/temp_extbnb_test/extbnb_tagger_results_stm_"    +std::to_string(random_number)+".txt", std::ios_base::app);
//		debug_file.open(path+"/temp_extbnb_test/extbnb_tagger_results_tgm_"    +std::to_string(random_number)+".txt", std::ios_base::app);
//		debug_file.open(path+"/temp_extbnb_test/extbnb_tagger_results_tgm_nof_"+std::to_string(random_number)+".txt", std::ios_base::app);

		debug_file << run_no << " " << subrun_no << " " << event_no << " " << tag_type << std::endl;
		debug_file.close();
		//std::cout << "run subrun event: " << run_no << " " << subrun_no << " " << event_no << ", tag_type = " << tag_type << std::endl;
	}

//	std::cout << run_no << " " << subrun_no << " " << event_no << " tag_type = " << tag_type << std::endl;
//	std::cout << "ending cosmic tagger ================================================================================" << std::endl;

	return std::make_tuple(tag_type, main_cluster, new_flash);
}

//Helper function that tells whether a set of extreme points is inside an x_bound
bool WCPPID::ToyFiducial::inside_x_region(std::vector<std::vector<WCP::WCPointCloud<double>::WCPoint>> extreme_points, double offset_x, double x_low, double x_high){
	bool inside = true;
	if(extreme_points[0][0].x-offset_x < x_low*units::cm || extreme_points[0][0].x-offset_x > x_high*units::cm
	|| extreme_points[1][0].x-offset_x < x_low*units::cm || extreme_points[1][0].x-offset_x > x_high*units::cm){
		inside = false;
	}
	return inside;
}

//Helper function that creates a vector of predicted PMT responses
std::vector<double> WCPPID::ToyFiducial::calculate_pred_pe(int run_no, double offset_x, WCP::Photon_Library *pl, WCPPID::PR3DCluster* main_cluster, std::vector<WCPPID::PR3DCluster*> additional_clusters, WCP::Opflash* flash, bool flag_data){

	std::vector<double> pred_pmt_light(32,0);
	double rel_light_yield_err = pl->rel_light_yield_err;
	double scaling_light_mag = pl->scaling_light_mag;
	std::map<int,int>* map_lib_pmt = &pl->map_lib_pmt;
	std::map<int,int>* map_pmt_lib = &pl->map_pmt_lib;
	std::vector<std::list<std::pair<int,float>>>* photon_library = &pl->library;

	double high_x_cut = 256 * units::cm;
	double high_x_cut_ext1 = + 1.2*units::cm;
	double high_x_cut_ext2 = - 2.0*units::cm;
	double low_x_cut = 0*units::cm;
	double low_x_cut_ext1 = - 2*units::cm;
	double low_x_cut_ext2 = + 4.0*units::cm;

	//Opflash *flash = bundle->get_flash();
	//double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
  	//PR3DCluster* main_cluster = bundle->get_main_cluster();

  	double first_pos_x = (*((main_cluster->get_time_cells_set_map().begin())->second.begin()))->get_sampling_points().front().x;
  	double last_pos_x = (*((main_cluster->get_time_cells_set_map().rbegin())->second.begin()))->get_sampling_points().front().x;

	bool flag_spec_end = false;
	
	// improve the position code ... 
	if (first_pos_x - offset_x <= low_x_cut + low_x_cut_ext1 &&
	    first_pos_x - offset_x > low_x_cut - 120*units::cm ){
	  
	  std::map<int,SMGCSet>& time_cells_set_map = main_cluster->get_time_cells_set_map();
	  int num_mcells_outside = 0;
	  int num_time_slices_outside = 0;

	  int num_mcells_def_outside = 0;
	  
	  double prev_pos_x= first_pos_x;
	  double current_pos_x = first_pos_x;
	  for (auto it3 = time_cells_set_map.begin(); it3 != time_cells_set_map.end(); it3++){
	    current_pos_x = (*(it3->second.begin()))->get_sampling_points().front().x;


	      // if (flash->get_flash_id()==35&&abs(main_cluster->get_cluster_id()-5)<=0)
	      // std::cout << num_time_slices_outside<< " " << num_mcells_outside << "  " << (first_pos_x - offset_x)/units::cm << " " <<
	      // 	(prev_pos_x-offset_x)/units::cm << " " << (current_pos_x-offset_x)/units::cm << std::endl;
	    
	    if (current_pos_x -offset_x > low_x_cut + low_x_cut_ext1 && current_pos_x - prev_pos_x > 0.75*units::cm)
	      break;
	    if (num_time_slices_outside > 60) break;
	    
	    if (current_pos_x -offset_x < low_x_cut + low_x_cut_ext1) num_mcells_def_outside += it3->second.size();

	    num_time_slices_outside += 1;
	    num_mcells_outside += it3->second.size();
	    prev_pos_x = current_pos_x;
	  }
	  if (num_time_slices_outside <=36 && num_mcells_outside < 0.05*main_cluster->get_num_mcells()){
	    first_pos_x = current_pos_x;
	    if (num_time_slices_outside > 10 && fabs(current_pos_x - prev_pos_x)<10*units::cm)
	      flag_spec_end = true;
	  }else if (num_time_slices_outside <=60 && num_mcells_outside < 0.06*main_cluster->get_num_mcells() && fabs(current_pos_x - prev_pos_x)>10*units::cm){
	    first_pos_x = current_pos_x;
	  }else if (num_time_slices_outside <=25 && num_mcells_outside < 0.12 * main_cluster->get_num_mcells() && fabs(current_pos_x - prev_pos_x)>20*units::cm){
	    first_pos_x = current_pos_x;
	  }

	  if (num_mcells_def_outside < 0.0015 * main_cluster->get_num_mcells()&&num_mcells_def_outside>0)
	    first_pos_x = offset_x;
	  
	  //  if (flash->get_flash_id()==61&&main_cluster->get_cluster_id()==1)
	  // std::cout << num_mcells_outside << " " << main_cluster->get_num_mcells() << "  A " << (first_pos_x - offset_x)/units::cm << " " <<
	  //  (prev_pos_x-offset_x)/units::cm << " " << (current_pos_x-offset_x)/units::cm << " " << num_mcells_def_outside << std::endl;
	  
	}
	if (last_pos_x - offset_x >= high_x_cut + high_x_cut_ext1 &&
	    last_pos_x - offset_x < high_x_cut + 120*units::cm){
	  std::map<int,SMGCSet>& time_cells_set_map = main_cluster->get_time_cells_set_map();
	  int num_mcells_outside = 0;
	  int num_time_slices_outside = 0;
	  int num_mcells_def_outside = 0;
	  double prev_pos_x= last_pos_x;
	  double current_pos_x = last_pos_x;

	  for (auto it3 = time_cells_set_map.rbegin(); it3 != time_cells_set_map.rend(); it3++){
	    current_pos_x = (*(it3->second.begin()))->get_sampling_points().front().x;
	    if (current_pos_x -offset_x<high_x_cut + high_x_cut_ext1 && fabs(current_pos_x - prev_pos_x) > 0.75*units::cm)
	      break;
	    if (num_time_slices_outside > 60) break;


	    if (current_pos_x -offset_x>high_x_cut + high_x_cut_ext1) num_mcells_def_outside +=it3->second.size();
	    
	    num_time_slices_outside += 1;
	    num_mcells_outside += it3->second.size();
	    prev_pos_x = current_pos_x;
	  }
	  if (num_time_slices_outside <=36 && num_mcells_outside < 0.05*main_cluster->get_num_mcells()){
	    last_pos_x = current_pos_x;
	    if (num_time_slices_outside > 10 && fabs(current_pos_x - prev_pos_x)<10*units::cm)
	      flag_spec_end = true;
	  }else if (num_time_slices_outside <=60 && num_mcells_outside < 0.06*main_cluster->get_num_mcells() && fabs(current_pos_x - prev_pos_x)>10*units::cm){
	    last_pos_x = current_pos_x;
	  }else if (num_time_slices_outside <=25 && num_mcells_outside < 0.12 * main_cluster->get_num_mcells() && fabs(current_pos_x - prev_pos_x)>20*units::cm){
	    last_pos_x = current_pos_x;
	  }
	  
	  if (num_mcells_def_outside < 0.0015 * main_cluster->get_num_mcells()&&num_mcells_def_outside>0)
	    last_pos_x = offset_x+high_x_cut;

	  // if (flash->get_flash_id()==19&&main_cluster->get_cluster_id()==19)
	  //   std::cout << flash->get_flash_id() << " "<< main_cluster->get_cluster_id() << " " << (first_pos_x-offset_x)/units::cm << " " << (last_pos_x-offset_x)/units::cm << " " << num_time_slices_outside << " " << num_mcells_outside << " " << main_cluster->get_num_mcells() << " " << fabs(current_pos_x - prev_pos_x)/units::cm << std::endl;
	  
	}
	
	// if (flash->get_flash_id()==16 && main_cluster->get_cluster_id()==13 )
	//   std::cout << flash->get_flash_id() << " "<< main_cluster->get_cluster_id() << " " << (first_pos_x-offset_x)/units::cm << " " << (last_pos_x-offset_x)/units::cm << std::endl;

	// if (flash->get_flash_id()==39)
	//   std::cout << flash->get_flash_id() << " " << main_cluster->get_cluster_id() << " " << offset_x/units::cm << " " << (first_pos_x-offset_x)/units::cm << " " << (last_pos_x-offset_x)/units::cm << " " << std::endl;
	
  	if (first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 -1.0*units::cm &&
  	    last_pos_x-offset_x > low_x_cut &&
  	    last_pos_x-offset_x < high_x_cut + high_x_cut_ext1 &&
  	    first_pos_x-offset_x < high_x_cut){
	  
	  //bundle->set_spec_end_flag(flag_spec_end);
	  
  	  //if (first_pos_x-offset_x <=low_x_cut + low_x_cut_ext2 && first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 - 1.0*units::cm ){
  	    // bundle->set_flag_close_to_PMT(true);
  	  //   bundle->set_flag_at_x_boundary(true);
  	  // }else if (first_pos_x-offset_x <= low_x_cut + low_x_cut_ext1 && first_pos_x-offset_x > low_x_cut + low_x_cut_ext1 -1.0*units::cm){
	    //bundle->set_flag_close_to_PMT(true);
  	    //bundle->set_flag_at_x_boundary(true);
	  //}
  	  //if (last_pos_x-offset_x >= high_x_cut + high_x_cut_ext2 && last_pos_x-offset_x < high_x_cut + high_x_cut_ext1){
  	    //bundle->set_flag_at_x_boundary(true);
  	  //}

  	  PR3DClusterSelection temp_clusters;
  	  temp_clusters.push_back(main_cluster);
  	  for (auto it3 = additional_clusters.begin(); it3!=additional_clusters.end(); it3++){
  	    temp_clusters.push_back(*it3);
  	  //  other_clusters->push_back(it3->first);
  	  }

  	  for (auto it3 = temp_clusters.begin(); it3!=temp_clusters.end(); it3++){
  	    SMGCSelection& mcells = (*it3)->get_mcells();
  	    bool flag_save = true;

	    if ((*it3) == main_cluster)  flag_save = false;
	    
  	    for (auto it4 = mcells.begin(); it4!=mcells.end(); it4++){
  	      SlimMergeGeomCell *mcell = (*it4);
  	      if (mcell->get_q()>0){
  		PointVector& pts = mcell->get_sampling_points();
  		if (pts.at(0).x-offset_x < low_x_cut+low_x_cut_ext1 ||
  		    pts.at(0).x-offset_x > high_x_cut+high_x_cut_ext1){
  		  flag_save = false;
  		  continue;
  		}
		
  		float charge = mcell->get_q()/pts.size();
  		Point p;
  		for (size_t i=0;i!=pts.size();i++){
  		  p.x = pts.at(i).x - offset_x;
  		  p.y = pts.at(i).y;
  		  p.z = pts.at(i).z;
		  
		  int voxel_id = convert_xyz_voxel_id(p);
		  std::list<std::pair<int,float>>& pmt_list = photon_library->at(voxel_id);
		  
		  for (auto it5 = pmt_list.begin(); it5!=pmt_list.end(); it5++){
		    pred_pmt_light[(*map_lib_pmt)[it5->first]] += charge * it5->second;
		  }
		}
  	      }
  	    }
  	    //if (flag_save)
  	    //  more_clusters->push_back(*it3);
	    
  	  } // loop over all clusters within this TPC object...

	  // veto PMT 18 ([17])
	  double norm_factor[32];
	  for (int i=0;i!=32;i++){
	    norm_factor[i] = 1;
	  }
	  if (flag_data){
	    if (run_no >= 12809)
	      norm_factor[17] = 0;
	  }
	  
  	  double sum1 = 0, sum2 = 0, max_pe = 0;
  	  for (size_t i=0;i!=32;i++){
  	    pred_pmt_light[i] *= scaling_light_mag * norm_factor[i];

  	    sum1 += flash->get_PE(i);
  	    sum2 += pred_pmt_light[i];
  	    if (pred_pmt_light[i] > max_pe)
  	      max_pe = pred_pmt_light[i];
  	  }
	  
	  // if (sum2 < sum1 * 3){ // three times allowrance ... 
	  //flag_good_bundle = true;
	  // }
  	}
	return pred_pmt_light;
}

//Helper Function
int WCPPID::ToyFiducial::convert_xyz_voxel_id(WCP::Point &p){
  // int voxel_x_id = std::round((p.x/units::cm+64.825-5.14667/2.)/5.14667);
  // int voxel_y_id = std::round((p.y/units::cm+193-5.14667/2.)/5.14667);
  // int voxel_z_id = std::round((p.z/units::cm+128.243-3.23122/2.)/3.23122);
  
  int voxel_x_id = std::round((p.x/units::cm+63.435-5.1096/2.)/5.1095);
  int voxel_y_id = std::round((p.y/units::cm+191.61-5.1096/2.)/5.1096);
  int voxel_z_id = std::round((p.z/units::cm+92.375-3.05437/2.)/3.05437);
  //fit
  // int voxel_y_id = std::round((p.y/units::cm+188.416-5.11698/2.)/5.11698);
  // int voxel_z_id = std::round((p.z/units::cm+88.3554-3.05768/2.)/3.05768);
  
  if (voxel_x_id<0) voxel_x_id=0;
  if (voxel_x_id>=75) voxel_x_id=74;
  if (voxel_y_id<0) voxel_y_id=0;
  if (voxel_y_id>=75) voxel_y_id=74;
  if (voxel_z_id<0) voxel_z_id=0;
  if (voxel_z_id>=400) voxel_z_id=399;

  int voxel_id = voxel_z_id*75*75 + voxel_y_id*75 + voxel_x_id;
  return voxel_id;
} 

//Helper function that returns the number of boundary contacts for a set of extreme points, drift offset, and tolerances
//checks whether a given cluster's extreme points are touching any detector boundaries for a given flash time / drift offest_x
//-2 means short track, -1 means outside boundary, 0 means inside boudnary, 1 means STM, 2 meanst TGM
int WCPPID::ToyFiducial::check_boundary(std::vector<std::vector<WCPointCloud<double>::WCPoint>> extreme_points, double offset_x, std::vector<double>* tol_vec){

	//Check whether the extreme points are contained within the boundary planes (allowing for tolerance), and if they are near a boundary.
	bool front_flag = false;
	bool back_flag = false;
	for(int i=0;i<int(extreme_points.size());i++){
		for(int ii=0;ii<int(extreme_points[i].size());ii++){
			WCP::Point p(extreme_points[i][ii].x,extreme_points[i][ii].y,extreme_points[i][ii].z);
			//Check whether the point is inside the extended fiducial volume
			std::vector<double> neg_tol_vec = {-1*tol_vec->at(0),-1*tol_vec->at(1),-1*tol_vec->at(2),-1*tol_vec->at(3), -1*tol_vec->at(4)};
			if(!inside_fiducial_volume(p,offset_x,tol_vec)){
				return -1;
			//Now that the point is known to be within the extended fiducial volume, check whether it is near any TGM boundaries, but only if it is near a PCA endpoint
			} else {
				if(!inside_fiducial_volume(p,offset_x,&neg_tol_vec)){
					if(i==0){front_flag = true;}
					if(i==1){back_flag = true;}
				}
			}
		}
	}
	return front_flag + back_flag;
}





//Enabled the fully contained check
bool WCPPID::ToyFiducial::check_fully_contained(WCPPID::PR3DCluster* main_cluster, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud){
  
  WCPPID::PR3DCluster *main_cluster1 = main_cluster;
  

  std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster1->get_extreme_wcps();

  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  Vector main_dir = main_cluster->get_PCA_axis(0);
  TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
  
  for (size_t i=0;i!=out_vec_wcps.size();i++){
    // check all the points ... 
    for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      if (!inside_fiducial_volume(p1,offset_x)){
        return false;
      }
    }


    Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
    TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
    dir *= (-1);

    // check U and V and W
    TVector3 dir_1(0,dir.Y(),dir.Z());
    double angle1 = dir_1.Angle(U_dir);
    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;
    
    double angle2 = dir_1.Angle(V_dir);
    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;
    
    double angle3 = dir_1.Angle(W_dir);
    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;
    
    // not added for now, need to check to add this one in when more events are available ...
    // XQ, 7/11/2018
    // double angle4 = fabs(3.1415926/2.-dir.Angle(drift_dir))/3.1415926*180.;


    //std::cout << "A: " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
      if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
        return false;
      }
    }
	    
    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
      if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	return false;
      }
    }
    
  }

  
  return true;
}


// check TGM code ...

bool WCPPID::ToyFiducial::check_tgm(WCPPID::PR3DCluster* main_cluster,Opflash* main_flash, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud){

  WCPPID::PR3DCluster *main_cluster1 = main_cluster; 
  std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster1->get_extreme_wcps();


  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  double length_limit = sqrt(pow(out_vec_wcps.at(0).at(0).x-out_vec_wcps.at(1).at(0).x,2)+
			     pow(out_vec_wcps.at(0).at(0).y-out_vec_wcps.at(1).at(0).y,2)+
			     pow(out_vec_wcps.at(0).at(0).z-out_vec_wcps.at(1).at(0).z,2));
  
  // std::cout << "Flash: " << flash->get_flash_id() << std::endl;

  // take a look at the first point ...
  for (size_t i=0;i!=out_vec_wcps.size();i++){
    bool flag_p1_inside = true;
    int p1_index = -1;
    for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      flag_p1_inside = flag_p1_inside && inside_fiducial_volume(p1,offset_x);
      if (!flag_p1_inside){
	p1_index = j;
	break;
      }
    }
    
    
    // loop through the remaining groups and check ...
    for (size_t k=i+1;k!=out_vec_wcps.size();k++){
      bool flag_p2_inside = true;
      int p2_index = -1;
      for(size_t j=0;j!=out_vec_wcps.at(k).size();j++){
	Point p2(out_vec_wcps.at(k).at(j).x,out_vec_wcps.at(k).at(j).y,out_vec_wcps.at(k).at(j).z);
	flag_p2_inside = flag_p2_inside && inside_fiducial_volume(p2,offset_x);
	if (!flag_p2_inside){
	  p2_index = j;
	  break;
	}
      }
      

      // if (main_cluster->get_cluster_id()==7){
      //  	std::cout << main_cluster->get_cluster_id() << " " << i << " " <<
      //  	  out_vec_wcps.at(i).at(0).x/units::cm << " " << out_vec_wcps.at(i).at(0).y/units::cm << " " << out_vec_wcps.at(i).at(0).z/units::cm << " " << 
      //  	  k << " " << out_vec_wcps.at(k).at(0).x/units::cm << " " << out_vec_wcps.at(k).at(0).y/units::cm << " " << out_vec_wcps.at(k).at(0).z/units::cm <<
      //  	  " " << p1_index << " " << p2_index << " " << flag_p1_inside << " " << flag_p2_inside << std::endl;
      // 	// for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
      // 	//   Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
      // 	//   std::cout << j << " A " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << std::endl;
      // 	// }
      // 	// for(size_t j=0;j!=out_vec_wcps.at(k).size();j++){
      // 	//   Point p2(out_vec_wcps.at(k).at(j).x,out_vec_wcps.at(k).at(j).y,out_vec_wcps.at(k).at(j).z);
      // 	//   std::cout << j << " B " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << std::endl;
      // 	// }
      // }
	
      if ((!flag_p1_inside) && (!flag_p2_inside)){
	// if not a neutrino candidate ... to be worked out ...
	 // Point p1(out_vec_wcps.at(i).at(p1_index).x,out_vec_wcps.at(i).at(p1_index).y,out_vec_wcps.at(i).at(p1_index).z);
	 // Point p2(out_vec_wcps.at(k).at(p2_index).x,out_vec_wcps.at(k).at(p2_index).y,out_vec_wcps.at(k).at(p2_index).z);
	//	std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << inside_fiducial_volume(p1,offset_x) << " A " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << inside_fiducial_volume(p2,offset_x) << " " << offset_x/units::cm << std::endl;


	// std::cout << main_cluster->get_cluster_id() << " " << (out_vec_wcps.at(i).at(p1_index).x-offset_x)/units::cm <<3=r-t=-[\] " " << out_vec_wcps.at(i).at(p1_index).y/units::cm << " " << out_vec_wcps.at(i).at(p1_index).z/units::cm << " " ;
	// std::cout << (out_vec_wcps.at(k).at(p2_index).x-offset_x)/units::cm << " " << out_vec_wcps.at(k).at(p2_index).y/units::cm << " " << out_vec_wcps.at(k).at(p2_index).z/units::cm << " " <<  flag_p1_inside << " " << flag_p2_inside << " " << out_vec_wcps.size() << std::endl;

	//std::cout << p1_index << " " << p2_index << std::endl;
	
	// check two points in between
	bool flag_check = false;
	for (int kk=0;kk!=3;kk++){
	  Point p3(out_vec_wcps.at(i).at(p1_index).x+ (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(p1_index).x)/4.*(kk+1),
		   out_vec_wcps.at(i).at(p1_index).y+ (out_vec_wcps.at(k).at(p2_index).y - out_vec_wcps.at(i).at(p1_index).y)/4.*(kk+1),
		   out_vec_wcps.at(i).at(p1_index).z+ (out_vec_wcps.at(k).at(p2_index).z - out_vec_wcps.at(i).at(p1_index).z)/4.*(kk+1));
	  flag_check = flag_check || inside_fiducial_volume(p3,offset_x);
	}
	
	// if (main_cluster->get_cluster_id()==10)
	//   std::cout << flag_check << " " << out_vec_wcps.size() << std::endl;
	
	if (flag_check){
	  if (main_flash->get_type()==2){
		    
	    double temp_length = sqrt(pow(out_vec_wcps.at(i).at(p1_index).x-out_vec_wcps.at(k).at(p2_index).x,2)+
				      pow(out_vec_wcps.at(i).at(p1_index).y-out_vec_wcps.at(k).at(p2_index).y,2)+
				      pow(out_vec_wcps.at(i).at(p1_index).z-out_vec_wcps.at(k).at(p2_index).z,2));

	    //std::cout << temp_length/units::cm << " " << length_limit/units::cm << std::endl;

	    if (i==0&&k==1){
	      if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud,true)) &&
		   temp_length > 0.45*length_limit)
		return true;
	    }else{
	      if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud)) &&
		   temp_length > 0.45*length_limit)
		return true;
	    }
	    
	  }else{
	    return true; // through going muon ...
	  }
	}else{
	  
	  if (out_vec_wcps.size()==2){
	    return true;
	  }else{
	    // if (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x))
	    //   return true;

	    bool flag_check_again = false;
	    for (int kkk = 0;kkk!=out_vec_wcps.size(); kkk++){
	      if (kkk==i || kkk==k) continue;
	      for (int kk=0;kk!=4;kk++){
	    	Point p3(out_vec_wcps.at(i).at(p1_index).x+ (out_vec_wcps.at(kkk).at(0).x - out_vec_wcps.at(i).at(p1_index).x)/4.*(kk+1),
	    		 out_vec_wcps.at(i).at(p1_index).y+ (out_vec_wcps.at(kkk).at(0).y - out_vec_wcps.at(i).at(p1_index).y)/4.*(kk+1),
	    		 out_vec_wcps.at(i).at(p1_index).z+ (out_vec_wcps.at(kkk).at(0).z - out_vec_wcps.at(i).at(p1_index).z)/4.*(kk+1));
	    	flag_check_again = flag_check_again || inside_fiducial_volume(p3,offset_x);
	      }
	      
	      for (int kk=0;kk!=3;kk++){
	    	Point p3(out_vec_wcps.at(kkk).at(0).x+ (out_vec_wcps.at(k).at(p2_index).x - out_vec_wcps.at(i).at(0).x)/4.*(kk+1),
	    		 out_vec_wcps.at(kkk).at(0).y+ (out_vec_wcps.at(k).at(p2_index).y - out_vec_wcps.at(i).at(0).y)/4.*(kk+1),
	    		 out_vec_wcps.at(kkk).at(0).z+ (out_vec_wcps.at(k).at(p2_index).z - out_vec_wcps.at(i).at(0).z)/4.*(kk+1));
	    	flag_check_again = flag_check_again || inside_fiducial_volume(p3,offset_x);
	      }
	    }
	    if (!flag_check_again){
	      //find the longest one ...
	      double temp_length = sqrt(pow(out_vec_wcps.at(i).at(p1_index).x-out_vec_wcps.at(k).at(p2_index).x,2)+
					pow(out_vec_wcps.at(i).at(p1_index).y-out_vec_wcps.at(k).at(p2_index).y,2)+
					pow(out_vec_wcps.at(i).at(p1_index).z-out_vec_wcps.at(k).at(p2_index).z,2));
	      
	      if (i==0&&k==1){
		if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud,true)) && temp_length > 0.45*length_limit)
		  return true;
	      }else{
		if ( (!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(p1_index),out_vec_wcps.at(k).at(p2_index),offset_x,ct_point_cloud)) && temp_length > 0.45*length_limit)
		  return true;
	      }
	    }	    
	  }
	}
      }else{
	Vector main_dir = main_cluster->get_PCA_axis(0);
	TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
	TVector3 dir_test(out_vec_wcps.at(i).at(0).x-out_vec_wcps.at(k).at(0).x,
			  out_vec_wcps.at(i).at(0).y-out_vec_wcps.at(k).at(0).y,
			  out_vec_wcps.at(i).at(0).z-out_vec_wcps.at(k).at(0).z);

	//	std::cout << main_cluster->get_cluster_id() << " " << fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.) << std::endl;

	// std::cout << main_cluster->get_cluster_id() << " " << (out_vec_wcps.at(i).at(0).x-offset_x)/units::cm << " " << out_vec_wcps.at(i).at(0).y/units::cm << " " << out_vec_wcps.at(i).at(0).z/units::cm << " " ;
	// std::cout << (out_vec_wcps.at(k).at(0).x-offset_x)/units::cm << " " << out_vec_wcps.at(k).at(0).y/units::cm << " " << out_vec_wcps.at(k).at(0).z/units::cm << " " <<  flag_p1_inside << " " << flag_p2_inside << " " << out_vec_wcps.size() << std::endl;
	
	if (fabs((3.1415926/2.-dir_test.Angle(dir_main))/3.1415926*180.)>75 || i==0 && k==1){
	  // check dead region ...
	  bool flag_p1_inside_p = flag_p1_inside;
	  if (flag_p1_inside_p){
	    Point p1(out_vec_wcps.at(i).at(0).x,out_vec_wcps.at(i).at(0).y,out_vec_wcps.at(i).at(0).z);
	    TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
	    dir *= (-1);

	    if (dir.Angle(dir_test) > 3.1415926*2./3.) continue;
	    
	    //	    std::cout << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;

	    // check U and V and W
	    TVector3 dir_1(0,dir.Y(),dir.Z());
	    double angle1 = dir_1.Angle(U_dir);
	    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
	    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

	    double angle2 = dir_1.Angle(V_dir);
	    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
	    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;

	    double angle3 = dir_1.Angle(W_dir);
	    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
	    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

	    // not added for now, need to check to add this one in when more events are available ...
	    // XQ, 7/11/2018
	    double angle4 = fabs(3.1415926/2.-dir.Angle(drift_dir))/3.1415926*180.;


	    //std::cout << "A: " << p1.x/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	      flag_p1_inside_p = flag_p1_inside_p && check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x);
	    }
	    
	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
	      flag_p1_inside_p= flag_p1_inside_p && check_dead_volume(p1,dir,1*units::cm,offset_x);
	  }
	  
	  bool flag_p2_inside_p = flag_p2_inside;
	  if (flag_p2_inside_p){
	    Point p2(out_vec_wcps.at(k).at(0).x,out_vec_wcps.at(k).at(0).y,out_vec_wcps.at(k).at(0).z);
	    TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
	    dir *= (-1);

	    if (dir.Angle(dir_test) < 3.1415926/3.) continue;
	    
	    //	    std::cout << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << dir.X() << " " << dir.Y() << dir.Z() << std::endl;
	    TVector3 dir_1(0,dir.Y(),dir.Z());
	    double angle1 = dir_1.Angle(U_dir);
	    TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
	    double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

	    double angle2 = dir_1.Angle(V_dir);
	    TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
	    double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;

	    double angle3 = dir_1.Angle(W_dir);
	    TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
	    double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

	    //	    std::cout << "B: " << p2.x/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " <<  angle1_1 << " " << angle2_1 << " " << angle3_1 << std::endl;

	    if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	      flag_p2_inside_p = flag_p2_inside_p && check_signal_processing(p2,dir,ct_point_cloud,1*units::cm,offset_x);
	    }
	    
	    if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 )
	      flag_p2_inside_p=flag_p2_inside_p && check_dead_volume(p2,dir,1*units::cm,offset_x);
	  }
	  
	  if ((!flag_p1_inside_p) && (!flag_p2_inside_p)){
	    if (main_flash->get_type()==2){
	      double temp_length = sqrt(pow(out_vec_wcps.at(i).at(0).x-out_vec_wcps.at(k).at(0).x,2)+
					pow(out_vec_wcps.at(i).at(0).y-out_vec_wcps.at(k).at(0).y,2)+
					pow(out_vec_wcps.at(i).at(0).z-out_vec_wcps.at(k).at(0).z,2));
	      if (i==0&&k==1){
		if ((!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(0),out_vec_wcps.at(k).at(0),offset_x,ct_point_cloud,true)) &&
		    temp_length > 0.45*length_limit
		    )
		  return true;
	      }else{
		if ((!check_neutrino_candidate(main_cluster1,out_vec_wcps.at(i).at(0),out_vec_wcps.at(k).at(0),offset_x,ct_point_cloud)) &&
		    temp_length > 0.45*length_limit
		    )
		  return true;
	      }
	    }else{
	      return true;
	    }
	  }
	  // check signal processing ...

	// {
	// 	if (flag_p1_inside)
	// 	  ;
	
	// 	if (flag_p2_inside)
	// 	  ;
	
	// }

	}
      }
    }
  }
  

  return false;

  // also check against the dead channel ...  
  // // check the fiducial volume ...
  // std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_main_axis_wcps();
  // Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
  // Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
  // //double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
  // bool flag_inside_p1 = inside_fiducial_volume(p1,offset_x);
  // bool flag_inside_p2 = inside_fiducial_volume(p2,offset_x);
  // //std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;

  
  
  // // check the dead region ... 
  // if (flag_inside_p1){
  //   // define a local direction ...
  //   TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
  //   dir *= (-1);
  //   flag_inside_p1=check_dead_volume(p1,dir,1*units::cm,offset_x);
  // }
  // if (flag_inside_p2){
  //   // define a  local direction ...
  //   TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
  //   dir *= (-1);
  //   flag_inside_p2=check_dead_volume(p2,dir,1*units::cm,offset_x);
  // }
  // return (!flag_inside_p1)&&(!flag_inside_p2);


  // bool flag_2nd = true;
       // {
	 
	 
       // 	 if ((!flag_inside_p1)&&(!flag_inside_p2)){
       // 	   event_type |= 1UL << 3; // through going muon ... 
       // 	   flag_2nd = false;
       // 	 }
	 
       // 	 if (flag_2nd && ((!flag_inside_p1)|| (!flag_inside_p2) )){
       // 	   // check the fiducial volume ...
       // 	   std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps = main_cluster->get_extreme_wcps();
	   
       // 	   Point p1(wcps.first.x,wcps.first.y,wcps.first.z);
       // 	   Point p2(wcps.second.x,wcps.second.y,wcps.second.z);
	   
       // 	   flag_inside_p1 = fid->inside_fiducial_volume(p1,offset_x);
       // 	   flag_inside_p2 = fid->inside_fiducial_volume(p2,offset_x);

       // 	   std::cout << main_cluster->get_cluster_id() << " " << (p1.x-offset_x)/units::cm << " " << p1.y/units::cm << " " << p1.z/units::cm << " " << (p2.x-offset_x)/units::cm << " " << p2.y/units::cm << " " << p2.z/units::cm << " " << fid->inside_fiducial_volume(p1,offset_x) << " " << fid->inside_fiducial_volume(p2,offset_x) << std::endl;
	   
       // 	   // check the dead region ...
       // 	   if (flag_inside_p1){
       // 	     // define a local direction ...
       // 	     TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
       // 	     dir *= (-1);
       // 	     flag_inside_p1=fid->check_dead_volume(p1,dir,1*units::cm,offset_x);
       // 	   }
       // 	   if (flag_inside_p2){
       // 	     // define a  local direction ...
       // 	     TVector3 dir = main_cluster->VHoughTrans(p2,30*units::cm);
       // 	     dir *= (-1);
       // 	     flag_inside_p2=fid->check_dead_volume(p2,dir,1*units::cm,offset_x);
       // 	   }
	   
       // 	   if ((!flag_inside_p1)&&(!flag_inside_p2)){
       // 	     event_type |= 1UL << 3; // through going muon ... 
       // 	   }
       // 	 }

  
  return false;
}



bool WCPPID::ToyFiducial::check_neutrino_candidate(WCPPID::PR3DCluster *main_cluster,WCPointCloud<double>::WCPoint& wcp1 ,WCPointCloud<double>::WCPoint& wcp2, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud,bool flag_2view_check){
  main_cluster->Create_graph(ct_point_cloud);
  
  main_cluster->dijkstra_shortest_paths(wcp1);
  main_cluster->cal_shortest_path(wcp2);

  std::list<WCPointCloud<double>::WCPoint>& path_wcps = main_cluster->get_path_wcps();
  
  
  
  PointVector path_wcps_vec;
  PointVector path_wcps_vec1;  
  double low_dis_limit = 0.5*units::cm;
  for (auto it = path_wcps.begin(); it!=path_wcps.end(); it++){
    if (path_wcps_vec.size()==0){
      Point p((*it).x,(*it).y,(*it).z);
      path_wcps_vec.push_back(p);
      path_wcps_vec1.push_back(p);
    }else{
      double dis = sqrt(pow((*it).x - path_wcps_vec.back().x,2)
			+pow((*it).y - path_wcps_vec.back().y,2)
			+pow((*it).z - path_wcps_vec.back().z,2));
      if (dis > low_dis_limit){
	Point p((*it).x,(*it).y,(*it).z);
	path_wcps_vec.push_back(p);
      }

      dis = sqrt(pow((*it).x - path_wcps_vec1.back().x,2)
		 +pow((*it).y - path_wcps_vec1.back().y,2)
		 +pow((*it).z - path_wcps_vec1.back().z,2));
      if (dis <= 2*low_dis_limit){
	Point p((*it).x,(*it).y,(*it).z);
	path_wcps_vec1.push_back(p);
      }else{
	int nseg = dis/2./low_dis_limit+1;
	for (int i=0;i!=nseg;i++){
	  Point temp_p;
	  temp_p.x = path_wcps_vec1.back().x + (i+1.)*((*it).x - path_wcps_vec1.back().x)/nseg;
	  temp_p.y = path_wcps_vec1.back().y + (i+1.)*((*it).y - path_wcps_vec1.back().y)/nseg;
	  temp_p.z = path_wcps_vec1.back().z + (i+1.)*((*it).z - path_wcps_vec1.back().z)/nseg;
	  path_wcps_vec1.push_back(temp_p);
	}
      }
    }
  }  

  
  //check whether path is good ... 
  {
    int num_nth = 0;
    double min_dis = 1e9;

    // bool flag_2view_check = true;
    
    // U and V induction view checks
    if (flag_2view_check){
      TVector3 drift_dir(1,0,0);
      // hard coded for U and V plane ... 
      TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
      TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
      TVector3 W_dir(0,1,0);

      TVector3 dir(wcp2.x-wcp1.x,wcp2.y-wcp1.y,wcp2.z-wcp1.z);
      
      TVector3 dir_1(0,dir.Y(),dir.Z());
      double angle1 = dir_1.Angle(U_dir);
      TVector3 tempV1(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle1),0);
      double angle1_1 = tempV1.Angle(drift_dir)/3.1415926*180.;

      double angle2 = dir_1.Angle(V_dir);
      TVector3 tempV2(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle2),0);
      double angle2_1 = tempV2.Angle(drift_dir)/3.1415926*180.;
      
      double angle3 = dir_1.Angle(W_dir);
      TVector3 tempV3(fabs(dir.X()), sqrt(dir.Y()*dir.Y()+dir.Z()*dir.Z())*sin(angle3),0);
      double angle3_1 = tempV3.Angle(drift_dir)/3.1415926*180.;

      double angle4 = fabs(3.1415926/2.-drift_dir.Angle(dir))/3.1415926*180.;
      
      
      if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5 || angle4 < 5.)){
	flag_2view_check = false;
      }
      
    }

    //int num_total_dead = 0;
    
    int num_bad = 0;
    
    for (int i=0;i!=path_wcps_vec1.size();i++){
      WCP::CTPointCloud<double> cloud_u = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,0);
      WCP::CTPointCloud<double> cloud_v = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,1);
      WCP::CTPointCloud<double> cloud_w = ct_point_cloud.get_closest_points(path_wcps_vec1.at(i),low_dis_limit*2,2);

      bool flag_reset = false;

      if (flag_2view_check){
	if (cloud_u.pts.size() >0 && cloud_v.pts.size() > 0 || // require two planes to be good ...
	    cloud_u.pts.size() >0 && cloud_w.pts.size() > 0 ||
	    cloud_v.pts.size() >0 && cloud_w.pts.size() > 0){
	  flag_reset = true;
	}else{
	  if (inside_dead_region(path_wcps_vec1.at(i)))
	    flag_reset = true;
	}
      }else{
	if (cloud_u.pts.size() >0  || // require one plane to be good ...
	    cloud_v.pts.size() >0  ||
	    cloud_w.pts.size() >0 ){
	  flag_reset = true;
	}else{
	  if (inside_dead_region(path_wcps_vec1.at(i)))
	    flag_reset = true;
	}
      }

      if (!ct_point_cloud.is_good_point(path_wcps_vec1.at(i)))
	num_bad ++;
      
      // std::cout << "O: " << path_wcps_vec1.at(i).x/units::cm << " " 
      // 		<< path_wcps_vec1.at(i).y/units::cm << " "
      // 		<< path_wcps_vec1.at(i).z/units::cm << " " << flag_reset << " " << cloud_u.pts.size() << " " << cloud_v.pts.size() << " "<< cloud_w.pts.size() << std::endl;
      
      if (flag_reset){
  	num_nth =0;
	min_dis = 1e9;
	num_bad = 0;
      }else{
	if (inside_fiducial_volume(path_wcps_vec1.at(i),offset_x)){
	  double dis1 = sqrt(pow(path_wcps_vec1.at(i).x-wcp1.x,2)+pow(path_wcps_vec1.at(i).y-wcp1.y,2)+pow(path_wcps_vec1.at(i).z-wcp1.z,2));
	  double dis2 = sqrt(pow(path_wcps_vec1.at(i).x-wcp2.x,2)+pow(path_wcps_vec1.at(i).y-wcp2.y,2)+pow(path_wcps_vec1.at(i).z-wcp2.z,2));
	  if (dis1 < min_dis) min_dis = dis1;
	  if (dis2 < min_dis) min_dis = dis2;
	  num_nth ++;
	  // num_total_dead ++;
	  
	  // if (main_cluster->get_cluster_id()==7)
	  //   std::cout << main_cluster->get_cluster_id() << " " << min_dis/units::cm << " " << flag_2view_check << " " << path_wcps_vec1.at(i).x/units::cm << " " 
	  // 	      << path_wcps_vec1.at(i).y/units::cm << " "
	  // 	      << path_wcps_vec1.at(i).z/units::cm << " "
	  // 	      << num_nth << " " << cloud_u.pts.size() << " " << cloud_v.pts.size() << " "<< cloud_w.pts.size() << " " << num_bad << std::endl;
	  
	}
	//std::cout << num_nth << std::endl;

	if (num_nth > 7 && min_dis < 25*units::cm && num_bad > 7) return true; // too big a gap ... 4 cm cut ...
	//if (num_nth > 7 && num_bad > 7) return true; // too big a gap ... 4 cm cut ...
      }
    }
  }
  
  
  // if (cluster_id == 13){
  //   std::cout << wcp1.x/units::cm << " " << wcp1.y/units::cm << " " << wcp1.z/units::cm << " " << wcp2.x/units::cm << " " << wcp2.y/units::cm << " " << wcp2.z/units::cm << std::endl;
  int count = 0;
  double max_angle=0;
  Point max_point(0,0,0);
  TVector3 drift_dir(1,0,0);
  for (size_t i=5;i+5<path_wcps_vec.size();i++){
    TVector3 dir1(path_wcps_vec.at(i).x - path_wcps_vec.at(i-5).x,
		  path_wcps_vec.at(i).y - path_wcps_vec.at(i-5).y,
		  path_wcps_vec.at(i).z - path_wcps_vec.at(i-5).z);
    TVector3 dir2(path_wcps_vec.at(i).x - path_wcps_vec.at(i+5).x,
		  path_wcps_vec.at(i).y - path_wcps_vec.at(i+5).y,
		  path_wcps_vec.at(i).z - path_wcps_vec.at(i+5).z);
    
    TVector3 dir3, dir4, dir5, dir6;
    {
      PointVector pts;
      double temp_x = 0;
      double temp_y = 0;
      double temp_z = 0;
      double temp_count = 0;
      for (size_t j=1;j!=15;j++){
       	if (i>=j){
	  Point pt(path_wcps_vec.at(i-j).x,path_wcps_vec.at(i-j).y,path_wcps_vec.at(i-j).z);
	  
	  if (j<=12&&j>2){
	    temp_x += pt.x;
	    temp_y += pt.y;
	    temp_z += pt.z;
	    temp_count ++;
	  }
	  pts.push_back(pt);
	}
	Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
	dir3 = main_cluster->calc_PCA_dir(pt,pts);
	dir5.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
		    temp_y/temp_count - path_wcps_vec.at(i).y,
		    temp_z/temp_count - path_wcps_vec.at(i).z);
	if (dir3.Angle(dir1)>3.1415926/2.)
	  dir3 *= -1;
      }
    }
    {
      PointVector pts;
      double temp_x = 0;
      double temp_y = 0;
      double temp_z = 0;
      double temp_count = 0;
      for (size_t j=1;j!=15;j++){
	if (i+j<path_wcps_vec.size()){
	  Point pt(path_wcps_vec.at(i+j).x,path_wcps_vec.at(i+j).y,path_wcps_vec.at(i+j).z);
	  if (j<=12&&j>2){
	    temp_x += pt.x;
	    temp_y += pt.y;
	    temp_z += pt.z;
	    temp_count ++;
	  }
	  pts.push_back(pt);
	}
      }
      Point pt(path_wcps_vec.at(i).x,path_wcps_vec.at(i).y,path_wcps_vec.at(i).z);
      dir4 = main_cluster->calc_PCA_dir(pt,pts);
      dir6.SetXYZ(temp_x/temp_count - path_wcps_vec.at(i).x,
      		  temp_y/temp_count - path_wcps_vec.at(i).y,
      		  temp_z/temp_count - path_wcps_vec.at(i).z);
      if (dir4.Angle(dir2)>3.1415926/2.)
	dir4 *= -1;
    }

    int cut1 = 0;
    if ((3.1415926 - dir1.Angle(dir2))/3.1415926*180.>25) cut1++;
    if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180.>25) cut1++;
    if ((3.1415926 - dir5.Angle(dir6))/3.1415926*180.>25) cut1++;
    int cut2 = 0;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. > 5) cut2++;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. > 5) cut2++;
    if (fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. > 5) cut2++;

    
    // if (main_cluster->get_cluster_id()==7)
    //   std::cout << i << " " << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << " " << (3.1415926 - dir1.Angle(dir2))/3.1415926*180. << " " << (3.1415926 - dir3.Angle(dir4))/3.1415926*180. << " " << (3.1415926 - dir5.Angle(dir6))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir1-dir2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir3-dir4))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(dir5-dir6))/3.1415926*180. << " " << cut1 << " " << cut2 << std::endl;
  
   
    
    
    if (cut1>=3 && cut2>=2){
      if ((3.1415926 - dir3.Angle(dir4))/3.1415926*180. > max_angle){
	max_angle = (3.1415926 - dir3.Angle(dir4))/3.1415926*180.;
	max_point = path_wcps_vec.at(i);
      }
      
      count ++;
      if (count >=3){
	TVector3 temp1(path_wcps_vec.at(i).x-wcp1.x,
		       path_wcps_vec.at(i).y-wcp1.y,
		       path_wcps_vec.at(i).z-wcp1.z);
	TVector3 temp2(path_wcps_vec.at(i).x-wcp2.x,
		       path_wcps_vec.at(i).y-wcp2.y,
		       path_wcps_vec.at(i).z-wcp2.z);

	// if (main_cluster->get_cluster_id()==7)
	//   std::cout << "A: " << (3.1415926-temp1.Angle(temp2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180.<< " " << temp1.Mag()/units::cm << " " << temp2.Mag()/units::cm << std::endl;


	if (((3.1415926-temp1.Angle(temp2))/3.1415926*180. >35 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5|| (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 ) ||
	    ((3.1415926-temp1.Angle(temp2))/3.1415926*180. >32 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5|| (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )&& temp1.Mag()>10*units::cm && temp2.Mag()>10*units::cm ||
	    ((3.1415926-temp1.Angle(temp2))/3.1415926*180. >25 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5 || (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )
	    && temp1.Mag()>15*units::cm && temp2.Mag()>15*units::cm){

	  // if (main_cluster->get_cluster_id()==7)
	  //   std::cout << "B: " <<  (!inside_fiducial_volume(max_point,offset_x)) << " " << inside_dead_region(max_point) << std::endl;
	  
	  if ((!inside_fiducial_volume(max_point,offset_x)) || // must be in fiducial
	      inside_dead_region(max_point)&&(3.1415926-temp1.Angle(temp2))/3.1415926*180<45 // not in dead_volume
	      ){ // should not too close to anode 
	  }else{
	    return true;
	  }
	}
      } else if (count>=1){
      	TVector3 temp1(path_wcps_vec.at(i).x-wcp1.x,
      		       path_wcps_vec.at(i).y-wcp1.y,
      		       path_wcps_vec.at(i).z-wcp1.z);
      	TVector3 temp2(path_wcps_vec.at(i).x-wcp2.x,
      		       path_wcps_vec.at(i).y-wcp2.y,
      		       path_wcps_vec.at(i).z-wcp2.z);

	//	std::cout << "B: " << path_wcps_vec.at(i).x/units::cm << " " << path_wcps_vec.at(i).y/units::cm << " " << path_wcps_vec.at(i).z/units::cm << " " << (3.1415926-temp1.Angle(temp2))/3.1415926*180. << " " << fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. << " " << temp1.Mag()/units::cm << " " << temp2.Mag()/units::cm << std::endl;
	
      	if (((3.1415926-temp1.Angle(temp2))/3.1415926*180. >35 && fabs(3.1415926/2.-drift_dir.Angle(temp1+temp2))/3.1415926*180. > 5.5 ||
	     (3.1415926-temp1.Angle(temp2))/3.1415926*180. >60 )
	    && temp1.Mag()>5*units::cm && temp2.Mag()>5*units::cm 
	    ){
	  //	  std::cout << "AAA" << std::endl; 

      	  if ((!inside_fiducial_volume(max_point,offset_x)) || // must be in fiducial
      	      inside_dead_region(max_point)&&(3.1415926-temp1.Angle(temp2))/3.1415926*180<45 // not in dead_volume
      	      ){ // should not too close to anode 
      	  }else{
      	    return true;
      	  }
      	}
      }
    }else{
      count = 0 ;
      max_angle = 0;
      max_point.x = 0;
      max_point.y = 0;
      max_point.z = 0;
    }
  }
  
  return false;
}
