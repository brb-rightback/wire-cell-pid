 
std::tuple<int, WCPPID::PR3DCluster*, WCP::Opflash*> WCPPID::ToyFiducial::glm_tagger(double eventTime, WCP::OpflashSelection& flashes, WCPPID::PR3DCluster* main_cluster, std::vector<WCPPID::PR3DCluster*> additional_clusters, WCP::Opflash* main_flash, std::tuple<int, double, double, int>& bundle_info, WCP::Photon_Library *pl, int time_offset, int nrebin, float unit_dis, WCP::ToyCTPointCloud& ct_point_cloud, int run_no, int subrun_no, int event_no, bool fully_contained, bool flag_match_data, bool flag_timestamp, bool debug_tagger){


	std::cout << "starting glm tagger ================================================================" << std::endl;

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
	double chi2_stm_tol = 1;
	double chi2_tgm_tol = 3;
	double bundle_ks_stm_tol =     0.04;					//This is the tolerance below which a flash is considered a very good match and should not be considered
	double bundle_ks_tgm_tol =     0.03;
	double bundle_ks_tgm_nof_tol = 0.03;
	double stm_flash_tol = 1.2*units::m;
	double tgm_flash_tol = 1.6*units::m;
	double stm_pe_frac_tol = 1.3;
	double tgm_pe_frac_tol = 1.5;
	std::vector<double> stm_tol_vec =     {0.0, 2.0, 2.0, 2.5, 2.0};	//x_ano, x_cat, ybot, ytop, z
	std::vector<double> tgm_tol_vec =     {2.2, 2.8, 2.8, 2.8, 2.8};
	std::vector<double> tgm_nof_tol_vec = {0.8, 1.2, 1.2, 1.2, 1.2};

/*
	if(!fully_contained){
		chi2_stm_tol = 1;
		chi2_tgm_tol = 1;
		bundle_ks_stm_tol =     0.04;
		bundle_ks_tgm_tol =     0.06;
		bundle_ks_tgm_nof_tol = 0.06;
		stm_flash_tol = 1.2*units::m;
		tgm_flash_tol = 1.2*units::m;
		stm_pe_frac_tol = 1.3;
		tgm_pe_frac_tol = 1.3;
		stm_tol_vec =     {0.0, 1.5, 1.5, 2.0, 1.5};
		tgm_tol_vec =     {1.6, 2.0, 2.0, 2.0, 2.0};
		tgm_nof_tol_vec = {0.6, 0.9, 0.9, 0.9, 0.9};
	}
*/
	if(!fully_contained){
		chi2_stm_tol = 1;
		chi2_tgm_tol = 1;
		bundle_ks_stm_tol =     0.04;
		bundle_ks_tgm_tol =     0.05;
		bundle_ks_tgm_nof_tol = 0.05;
		stm_flash_tol = 1.0*units::m;
		tgm_flash_tol = 1.0*units::m;
		stm_pe_frac_tol = 1.3;
		tgm_pe_frac_tol = 1.3;
		stm_tol_vec =     {0.0, 1.5, 1.5, 2.0, 1.5};
		tgm_tol_vec =     {2.0, 2.5, 2.5, 2.5, 2.5};
		tgm_nof_tol_vec = {0.6, 0.9, 0.9, 0.9, 0.9};
	}

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

	//Check whether the track is at least 10 cm long
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
			std::vector<double> pred_pmt_light = calculate_pred_pe(eventTime, run_no, offset_x, pl, main_cluster, additional_clusters, flash, flag_match_data, flag_timestamp);

			//Compute the PE centroids for the flash PE and predicted PE
			double pred_pe_tot = 0;
			double flash_pe_tot = 0;
			double pred_pe_z_centroid = 0;
			double flash_pe_z_centroid = 0;
			double newf_chi2 = 0;
			int    newf_ndf = 0;
			for(int i=0;i<int(pred_pmt_light.size());i++){
				double flash_pe = flash->get_PE(i);
				double flash_pe_err = flash->get_PE_err(i);
				if(flash_pe_err != 0 && (flash_pe != 0 || pred_pmt_light[i] != 0)){
					newf_chi2 += pow((flash_pe-pred_pmt_light[i])/flash_pe_err,2);
					newf_ndf++;
				}
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
			bool cathode_side, reasonable_pe, nearby_flash, better_chi2;
			double ks_frac, chi2_frac;

			//TGM-with-flash case
			ks_frac = bundle_ks / bundle_ks_tgm_tol;
			chi2_frac = (newf_chi2/newf_ndf)/(bundle_chi2/bundle_ndf);
			better_chi2 = chi2_frac < chi2_tgm_tol;
			reasonable_pe = flash_pe_tot < pred_pe_tot*tgm_pe_frac_tol;
			nearby_flash = flash_pe_z_centroid > pred_pe_z_centroid-tgm_flash_tol && flash_pe_z_centroid < pred_pe_z_centroid+tgm_flash_tol;
			if(better_chi2 && reasonable_pe && nearby_flash){
				std::vector<double> temp_tgm_tol_vec = tgm_tol_vec;
				if(ks_frac > 2){for(int i=0;i<5;i++){temp_tgm_tol_vec[i] *= 1.3;}}
				if(ks_frac > 1){
					boundary_num_tgm = check_boundary(extreme_points, offset_x, &temp_tgm_tol_vec);
				}
			}

			//STM-with-flash case
			cathode_side = inside_x_region(extreme_points,offset_x,128-stm_tol_vec[0],256+stm_tol_vec[1]);
			ks_frac = bundle_ks / bundle_ks_stm_tol;
			better_chi2 = chi2_frac < chi2_stm_tol;
			reasonable_pe = flash_pe_tot < pred_pe_tot*stm_pe_frac_tol;
			nearby_flash = flash_pe_z_centroid > pred_pe_z_centroid-stm_flash_tol && flash_pe_z_centroid < pred_pe_z_centroid+stm_flash_tol;
			if(better_chi2 && cathode_side && reasonable_pe && nearby_flash){
				std::vector<double> temp_stm_tol_vec = stm_tol_vec;
				if(ks_frac > 3){for(int i=0;i<5;i++){temp_stm_tol_vec[i] *= 1.3;}}
				if(ks_frac > 2){for(int i=0;i<5;i++){temp_stm_tol_vec[i] *= 1.3;}}
				if(ks_frac > 1){
					boundary_num_stm = check_boundary(extreme_points, offset_x, &temp_stm_tol_vec);
				}
			}

			//TGM gets priority over STM
			if(boundary_num_tgm==2)					{boundary_num =  2;}
			else if(boundary_num_stm==1)				{boundary_num =  1;}
			else if(boundary_num_tgm<0 || boundary_num_stm<0)	{boundary_num = -1;}

			//Store a STM / TGM tag if appropriate
			bool tgm = check_tgm(main_cluster,flash,offset_x,ct_point_cloud);
			if(boundary_num==1 || boundary_num==2 && tgm){
				n_boundary_list.push_back(boundary_num);
				candidate_flash_list.push_back(flash);
			}
		}

		//Check whether a TGM w/o flash hypothesis works
		//Skip if the existing neutrino match is very good
		double ks_frac = bundle_ks / bundle_ks_tgm_nof_tol;
		if(ks_frac > 2){for(int i=0;i<5;i++){tgm_nof_tol_vec[i] *= 1.3;}}
		if(ks_frac > 1){
			double offset_x = (main_flash->get_time() - time_offset)*2./nrebin*time_slice_width;
			double tx     = tgm_nof_tol_vec[0];
			double ty_bot = tgm_nof_tol_vec[2];
			double ty_top = tgm_nof_tol_vec[3];
			double tz     = tgm_nof_tol_vec[4];
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

/*
	debug_tagger = false;
	//Write to file for debugging purposes
	if(debug_tagger){
		for(int n_match = 0; n_match<int(n_boundary_list.size());n_match++){
			if(n_match < int(candidate_flash_list.size())){
				std::cout << "match type = " << n_boundary_list[n_match] << ", flash # = " << candidate_flash_list[n_match]->get_flash_id() << std::endl;		
			} else {
				std::cout << "match type = " << n_boundary_list[n_match] << std::endl;
			}
		}
		std::cout << run_no << " " << subrun_no << " " << event_no << " tag_type = " << tag_type << std::endl;
		write_debug(run_no,subrun_no,event_no,tag_type);
	}
*/
	std::cout << run_no << " " << subrun_no << " " << event_no << " tag_type = " << tag_type << ", ks = " << bundle_ks << std::endl;
	std::cout << "ending glm tagger ==================================================================" << std::endl;
	return std::make_tuple(tag_type, main_cluster, new_flash);
}

//Write to file with this function
void WCPPID::ToyFiducial::write_debug(int run_no, int subrun_no, int event_no, int tag_type){

	srand(time(NULL));
	int random_number = std::rand();
	std::ofstream debug_file;
	std::string path = "/uboone/data/users/lcoopert/cosmic_tagger/data";
	debug_file.open(path+"/temp_numu_cc_all/numu_cc_tagger_results_"+std::to_string(random_number)+".txt", std::ios_base::app);
//	debug_file.open(path+"/temp_extbnb_all/extbnb_tagger_results_"+std::to_string(random_number)+".txt", std::ios_base::app);

//	debug_file.open(path+"/temp_extbnb_test/extbnb_tagger_results_stm_"    +std::to_string(random_number)+".txt", std::ios_base::app);
//	debug_file.open(path+"/temp_extbnb_test/extbnb_tagger_results_tgm_"    +std::to_string(random_number)+".txt", std::ios_base::app);
//	debug_file.open(path+"/temp_extbnb_test/extbnb_tagger_results_tgm_nof_"+std::to_string(random_number)+".txt", std::ios_base::app);

	debug_file << run_no << " " << subrun_no << " " << event_no << " " << tag_type << std::endl;
	debug_file.close();

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
std::vector<double> WCPPID::ToyFiducial::calculate_pred_pe(double eventTime, int run_no, double offset_x, WCP::Photon_Library *pl, WCPPID::PR3DCluster* main_cluster, std::vector<WCPPID::PR3DCluster*> additional_clusters, WCP::Opflash* flash, bool flag_match_data, bool flag_timestamp){

  std::cout << "ZXin_3: " << eventTime << " " << flag_timestamp << std::endl;
  
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
	  if (flag_match_data){
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
			std::vector<double>* neg_tol_vec = NULL;
			if(tol_vec!=NULL){
				std::vector<double> temp_neg_tol_vec = {-1*tol_vec->at(0),-1*tol_vec->at(1),-1*tol_vec->at(2),-1*tol_vec->at(3), -1*tol_vec->at(4)};
				neg_tol_vec = &temp_neg_tol_vec;
			}
			if(!inside_fiducial_volume(p,offset_x,tol_vec)){
				return -1;
			//Now that the point is known to be within the extended fiducial volume, check whether it is near any TGM boundaries, but only if it is near a PCA endpoint
			} else {
				if(!inside_fiducial_volume(p,offset_x,neg_tol_vec)){
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


//edited version of check_stm that ignores tgms
bool WCPPID::ToyFiducial::check_stm_only(WCPPID::PR3DCluster* main_cluster, std::vector<WCPPID::PR3DCluster*>& additional_clusters, double offset_x, double flash_time, WCP::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, int& event_type){

  //  check_full_detector_dead();

  


  //std::cout << flag_other_clusters << std::endl;
  
  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  Vector main_dir = main_cluster->get_PCA_axis(0);
  TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
  
  std::vector<WCPointCloud<double>::WCPoint> candidate_exit_wcps;
  std::set<int> temp_set;
  std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps;
  
  // first round check
  if (main_cluster->get_point_cloud_steiner()->get_cloud().pts.size()==0)
    return false;

  {
    // get two extreme points ...
    wcps = main_cluster->get_two_boundary_wcps(2,true);
    // figure out the end points ...
    std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster->get_extreme_wcps();
    
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.first);
      out_vec_wcps.push_back(temp_wcps);
    }
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.second);
      out_vec_wcps.push_back(temp_wcps);
    }
    
    // boundary check
    for (size_t i=0;i!=out_vec_wcps.size();i++){
      bool flag_save = false;
      // check all the points ... 
      for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){

	Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
	//std::cout << p1 << " " << inside_fiducial_volume(p1,offset_x) << std::endl;
	if (!inside_fiducial_volume(p1,offset_x)){
	  candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  flag_save = true;
	  break;
	}
      }
      
      if (!flag_save){
	// check direction
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
	
	
	if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	  if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
	    flag_save = true;
	    candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  }
	}
	
	if (!flag_save){
	  if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
	    if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	      flag_save = true;
	      candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	    }
	  }
	}
      }
    }
    
    //std::cout << wcps.first.x << " " <<wcps.first.y << " " << wcps.first.z << std::endl;
    //std::cout << wcps.second.x << " " <<wcps.second.y << " " << wcps.second.z << std::endl;
    
    for (size_t i=0; i!=candidate_exit_wcps.size(); i++){
      double dis1 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.first.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.first.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.first.z,2));
      double dis2 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.second.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.second.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.second.z,2));
      
      //      std::cout << candidate_exit_wcps.at(i).x << " " << candidate_exit_wcps.at(i).y << " " << candidate_exit_wcps.at(i).z << " " << dis1 << " " << dis2 << std::endl;
      
      // essentially one of the extreme points ...
      if (dis1 < dis2){
	if (dis1 < 1.0*units::cm)  temp_set.insert(0);
      }else{
	if (dis2 < 1.0*units::cm)  temp_set.insert(1);
      }
    }

    //    std::cout << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;
    
    // protection against two end point situation
    if (temp_set.size()==2){
      WCP::Point tp1(wcps.first.x,wcps.first.y,wcps.first.z);
      WCP::Point tp2(wcps.second.x,wcps.second.y,wcps.second.z);
      
      temp_set.clear();
      
      if ((!inside_fiducial_volume(tp1,offset_x))) temp_set.insert(0);
      if ((!inside_fiducial_volume(tp2,offset_x))) temp_set.insert(1);
      if (temp_set.size()==0){
	temp_set.insert(0);
	temp_set.insert(1);
      }
    }
  }

  if (temp_set.size()==0){
    candidate_exit_wcps.clear();
    // get two extreme points ...
    wcps = main_cluster->get_two_boundary_wcps(2);
    // figure out the end points ...
    std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster->get_extreme_wcps();
    
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.first);
      out_vec_wcps.push_back(temp_wcps);
    }
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.second);
      out_vec_wcps.push_back(temp_wcps);
    }
    
    // boundary check
    for (size_t i=0;i!=out_vec_wcps.size();i++){
      bool flag_save = false;
      // check all the points ... 
      for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
	Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
	if (!inside_fiducial_volume(p1,offset_x)){
	  candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  flag_save = true;
	  break;
	}
      }
      
      if (!flag_save){
	// check direction
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
	
	
	if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	  if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
	    flag_save = true;
	    candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  }
	}

	//	std::cout << fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.) << " " << check_dead_volume(p1,dir,1*units::cm,offset_x) << std::endl;
	if (!flag_save){
	  if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
	    if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	      flag_save = true;
	      candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	    }
	  }
	}
      }
    }
    
    //std::cout << wcps.first.x << " " <<wcps.first.y << " " << wcps.first.z << std::endl;
    //std::cout << wcps.second.x << " " <<wcps.second.y << " " << wcps.second.z << std::endl;
    
    for (size_t i=0; i!=candidate_exit_wcps.size(); i++){
      double dis1 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.first.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.first.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.first.z,2));
      double dis2 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.second.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.second.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.second.z,2));
      
      //std::cout << candidate_exit_wcps.at(i).x << " " << candidate_exit_wcps.at(i).y << " " << candidate_exit_wcps.at(i).z << " " << dis1 << " " << dis2 << std::endl;
      
      // essentially one of the extreme points ...
      if (dis1 < dis2){
	if (dis1 < 1.0*units::cm)  temp_set.insert(0);
      }else{
	if (dis2 < 1.0*units::cm)  temp_set.insert(1);
      }
    }

    //std::cout << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;
    
    
    // protection against two end point situation
    if (temp_set.size()==2){
      WCP::Point tp1(wcps.first.x,wcps.first.y,wcps.first.z);
      WCP::Point tp2(wcps.second.x,wcps.second.y,wcps.second.z);
      
      temp_set.clear();
      
      if ((!inside_fiducial_volume(tp1,offset_x))) temp_set.insert(0);
      if ((!inside_fiducial_volume(tp2,offset_x))) temp_set.insert(1);
      if (temp_set.size()==0){
	temp_set.insert(0);
	temp_set.insert(1);
      }
    }
  }
  

  
  // fully contained, so not a STM
  if (candidate_exit_wcps.size()==0) {
    std::cout << "Mid Point: A" << std::endl;
    return false;
  }
  
  std::cout << "end_point: " << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;

  // It is possible that we have two points out of fiducial
  // Michel electron is outside the boundary ...

  
  
  // Crawl backward according to the graph ??? ...
  WCPointCloud<double>::WCPoint first_wcp;
  WCPointCloud<double>::WCPoint last_wcp;
  bool flag_double_end = false;

  if (temp_set.size()!=0){
    if (*temp_set.begin()==0){
      first_wcp = wcps.first;
      last_wcp = wcps.second;
    }else{
      first_wcp = wcps.second;
      last_wcp = wcps.first;
    }
    if (temp_set.size()==2) flag_double_end = true;
    
  }else{
    if (candidate_exit_wcps.size()==1){
      first_wcp = candidate_exit_wcps.at(0);
      TVector3 dir1(wcps.first.x - candidate_exit_wcps.at(0).x,
		    wcps.first.y - candidate_exit_wcps.at(0).y,
		    wcps.first.z - candidate_exit_wcps.at(0).z);
      TVector3 dir2(wcps.second.x - candidate_exit_wcps.at(0).x,
		    wcps.second.y - candidate_exit_wcps.at(0).y,
		    wcps.second.z - candidate_exit_wcps.at(0).z);
      double dis1 = dir1.Mag();
      double dis2 = dir2.Mag();

      if (dir1.Angle(dir2) > 120/180.*3.1415926 && dis1 > 20*units::cm &&
      	  dis2 > 20*units::cm){
	std::cout << "Mid Point: B" << std::endl;
	return false;
      }else{
	if (dis1 < dis2){
	  last_wcp = wcps.second;	   
	}else{
	  last_wcp = wcps.first; 
	}
      }
      
    }else{
      std::cout << "Mid Point: C" << std::endl;
      return false;
    }
  }

 
  bool flag_other_clusters = check_other_clusters(main_cluster, additional_clusters);
  //  std::cout << "haha " << flag_other_clusters << std::endl;
  
  // forward check ...
  {
    if (flag_double_end) std::cout << "Forward check! " << std::endl;
    // regular crawling ...
    main_cluster->do_rough_path(first_wcp, last_wcp);
    // std::cout << "haha" << std::endl;
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    //std::cout << "haha" << std::endl;
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
    //std::cout << "test" << std::endl; 
    if (main_cluster->get_fine_tracking_path().size()<=3) return false;
    
    //std::cout << "haha " << std::endl;
    Point mid_p = main_cluster->adjust_rough_path(); 
    //std::cout << "haha" << std::endl;
    // fitting trajectory and dQ/dx...
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);

    //    std::cout << "test" << std::endl; 
    // check both end points for TGM ...
    WCP::PointVector& pts = main_cluster->get_fine_tracking_path();
    std::vector<double>& dQ = main_cluster->get_dQ();
    std::vector<double>& dx = main_cluster->get_dx();
    
    int kink_num = find_first_kink(main_cluster);
    
    double left_L = 0; 
    double left_Q = 0;
    double exit_L = 0; 
    double exit_Q = 0;
    for (size_t i=0;i!=kink_num;i++){
      exit_L += dx.at(i);
      exit_Q += dQ.at(i);
    }
    for (size_t i = kink_num; i!=dx.size(); i++){
      left_L += dx.at(i);
      left_Q += dQ.at(i);
    }
    
     std::cout << "Left: " << exit_L/units::cm << " " << left_L/units::cm << " " << (left_Q/(left_L/units::cm+1e-9))/50e3 << " " << (exit_Q/(exit_L/units::cm+1e-9)/50e3) << std::endl;

     //std::cout << pts.front() << " " << pts.back() << " " << (!inside_fiducial_volume(pts.front(),offset_x)) << " " << (!inside_fiducial_volume(pts.back(),offset_x)) << std::endl;
    if ( (!inside_fiducial_volume(pts.front(),offset_x)) && (!inside_fiducial_volume(pts.back(),offset_x))){
      // std::cout << exit_L_TGM/units::cm << " " << left_L_TGM/units::cm << std::endl;
      bool flag_TGM_anode = false;
      
      if ((pts.back().x < 2*units::cm || pts.front().x < 2*units::cm) && kink_num >=0 && kink_num < pts.size()){ // at Anode ...
	if (pts.at(kink_num).x < 6*units::cm){
	  TVector3 v10(pts.back().x-pts.at(kink_num).x,pts.back().y-pts.at(kink_num).y,pts.back().z-pts.at(kink_num).z);
	  TVector3 v20(pts.front().x-pts.at(kink_num).x,pts.front().y-pts.at(kink_num).y,pts.front().z-pts.at(kink_num).z);
	  if (fabs(v10.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v10.Mag()>15*units::cm || fabs(v20.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v20.Mag()>15*units::cm)
	      flag_TGM_anode = true;
	  //std::cout << v10.Angle(drift_dir)/3.1415926*180. << " " << v20.Angle(drift_dir)/3.1415926*180. << std::endl;
	}
      }
      if ((exit_L < 3*units::cm || left_L < 3*units::cm || flag_TGM_anode)){
	std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
//	event_type |= 1UL << 3;
//	return true;
      }
    }else if ((!inside_fiducial_volume(pts.front(),offset_x)) && left_L < 3*units::cm){
      // check dead volume and SP ...
      Point p1 = pts.back();
      TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
      dir *= (-1);
      if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	if (exit_L < 3*units::cm || left_L < 3*units::cm){
	  std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
//	  event_type |= 1UL << 3;
//	  return true;
	}
      }
      //std::cout << "ABC: " << dir.X() << " " << dir.Y() << " " << dir.Z() << " " << check_dead_volume(p1,dir,1*units::cm,offset_x) << std::endl;
    }
  
    if (left_L > 40*units::cm || left_L > 7.5*units::cm && (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0){
      if (!flag_double_end){
	std::cout << "Mid Point A " << inside_dead_region(mid_p) << " " << mid_p << " " << left_L  << " " << (left_Q/(left_L/units::cm+1e-9)/50e3) << std::endl;
	return false;
      }
    }else{
      bool flag_fix_end = false;
      if (exit_L < 35*units::cm || (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0&& left_L > 2*units::cm) flag_fix_end = true;
      
      if (left_L < 8*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3)< 1.5 ||
	  left_L < 6*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.7 ||
	  left_L < 5*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.8 ||
	  left_L < 3*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.9){
	left_L = 0;
	kink_num = dQ.size();
	exit_L = 40*units::cm;
	flag_fix_end = false;
      }
    
     
      bool flag_pass = false;

      // std::cout << flag_other_clusters << std::endl;
      
      if (!flag_other_clusters){
	if (left_L < 40*units::cm) {
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
	
	if (left_L < 20*units::cm){
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
	    }
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm , 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm , 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
      }else{
	if (flag_fix_end)
	  flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm, true);
	else
	  flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm, 0., 35*units::cm, true);
      }
      
      if (flag_pass) {
	main_cluster->clear_fit_tracks();
	main_cluster->search_other_tracks(ct_point_cloud, global_wc_map, flash_time*units::microsecond);

	if (check_other_tracks(main_cluster, offset_x)){
	  std::cout << "Mid Point Tracks" << std::endl;
	  return false;
	}
	
	//	std::cout << main_cluster->get_fit_tracks().size() << std::endl;
	
	if (!detect_proton(main_cluster, kink_num)) return true;
      }
    }
  }
  
  // backward check ...
  if (flag_double_end){
    std::cout << "Backward check! " << std::endl;
    // regular crawling ...
    main_cluster->do_rough_path( last_wcp, first_wcp);
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
    Point mid_p = main_cluster->adjust_rough_path(); 
    // fitting trajectory and dQ/dx...
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);

    // check both end points for TGM ...
    WCP::PointVector& pts = main_cluster->get_fine_tracking_path();
    std::vector<double>& dQ = main_cluster->get_dQ();
    std::vector<double>& dx = main_cluster->get_dx();
    
    int kink_num = find_first_kink(main_cluster);
    
    double left_L = 0;
    double left_Q = 0;
    double exit_L = 0;
    double exit_Q = 0;
    for (size_t i=0;i!=kink_num;i++){
      exit_L += dx.at(i);
      exit_Q += dQ.at(i);
    }
    for (size_t i = kink_num; i!=dx.size(); i++){
      left_L += dx.at(i);
      left_Q += dQ.at(i);
    }
    
    std::cout << "Left: " << exit_L/units::cm << " " << left_L/units::cm << " " << (left_Q/(left_L/units::cm+1e-9))/50e3 << " " << (exit_Q/(exit_L/units::cm+1e-9)/50e3) << std::endl;
    
    if ( (!inside_fiducial_volume(pts.front(),offset_x)) && (!inside_fiducial_volume(pts.back(),offset_x))){
      bool flag_TGM_anode = false;
      
      if ((pts.back().x < 2*units::cm || pts.front().x < 2*units::cm) && kink_num >=0 && kink_num < pts.size()){ // at Anode ...
	if (pts.at(kink_num).x < 6*units::cm){
	  TVector3 v10(pts.back().x-pts.at(kink_num).x,pts.back().y-pts.at(kink_num).y,pts.back().z-pts.at(kink_num).z);
	  TVector3 v20(pts.front().x-pts.at(kink_num).x,pts.front().y-pts.at(kink_num).y,pts.front().z-pts.at(kink_num).z);
	  if (fabs(v10.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v10.Mag()>15*units::cm || fabs(v20.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v20.Mag()>15*units::cm)
	      flag_TGM_anode = true;
	  //std::cout << v10.Angle(drift_dir)/3.1415926*180. << " " << v20.Angle(drift_dir)/3.1415926*180. << std::endl;
	}
      }
      if ((exit_L < 3*units::cm || left_L < 3*units::cm) || flag_TGM_anode){
	std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
//	event_type |= 1UL << 3;
//	return true;
      }
    }
  
    if (left_L > 40*units::cm || left_L > 7.5*units::cm && (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0){
      std::cout << "Mid Point A " << inside_dead_region(mid_p) << " " << mid_p << " " << left_L  << " " << (left_Q/(left_L/units::cm+1e-9)/50e3) << std::endl;
      return false;
    }else{
      bool flag_fix_end = false;
      if (exit_L < 35*units::cm || (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0 && left_L > 2*units::cm) flag_fix_end = true;
      
      if (left_L < 8*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3)< 1.5 ||
	  left_L < 6*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.7 ||
	  left_L < 3*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.9){
	left_L = 0;
	kink_num = dQ.size();
	exit_L = 40*units::cm;
	flag_fix_end = false;
      }

    
    
      bool flag_pass = false;
      if (!flag_other_clusters){
	if (left_L < 40*units::cm) {
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
	
	if (left_L < 20*units::cm){
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
	    }
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm , 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm , 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
      }else{
	if (flag_fix_end)
	  flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm, true);
	else
	  flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm, 0., 35*units::cm, true);
      }

      if (flag_pass) {
	main_cluster->clear_fit_tracks();
	main_cluster->search_other_tracks(ct_point_cloud, global_wc_map, flash_time*units::microsecond);
	if (check_other_tracks(main_cluster, offset_x)){
	  std::cout << "Mid Point Tracks" << std::endl;
	  return false;
	}
	
	if (!detect_proton(main_cluster, kink_num)) return true;
      }
    }
  }
  
  // check 5512-209-10491
  // if (inside_dead_region(mid_p)) return true;
  //  

  // end check ...
  std::cout << "Mid Point " << std::endl;  
  return false;
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



//////////////////////////////////////////////////// M2
//////////////////////////////////////////////////// M2

bool WCPPID::ToyFiducial::M2_distance_cut(WCPPID::PR3DCluster* main_cluster, double user_offset_dist, double user_offset_frac)
{
  bool result = false;

  //TVector3 vp0(p0.x/units::cm, p0.y/units::cm, p0.z/units::cm);
  //TVector3 vp1(p1.x/units::cm, p1.y/units::cm, p1.z/units::cm);
  //TVector3 vDir = vp1 - vp0;
  //vDir = 1./vDir.Mag() * vDir;
  
  //Vector main_dir = main_cluster->get_PCA_axis(0);
  //TVector3 vDir(main_dir.x,main_dir.y,main_dir.z);
  //int n_total = 0;
  //int n_offset = 0;
  
  return result;
}



double WCPPID::ToyFiducial::M2_offset_YX_x(WCP::Point& p)
{
  double result = 0;
    
  double eff_x = p.x;
  double eff_y = p.y;
  double eff_z = p.z;

  /////////////////////////

  int index_z = 100;

  if( eff_z/units::cm<0 ) {
    index_z = 1;
  }
  else if( eff_z/units::cm<=1037 ) {
    index_z = (int)(eff_z/units::cm/100)+1;
    if( index_z>10 ) index_z = 10;
  }
  else {
    index_z = 10;
  }

  /*
    ---------\   TOP-1
              \
               \ TOP-2
               |
               |
               |
               |
               | BOT-2
               /
              /
    ---------/   BOT-1
  */

  //std::cout<<" xp test p "<<eff_x/units::cm<<" "<<eff_y/units::cm<<" "<<eff_z/units::cm<<std::endl;
  //std::cout<<" xp test   "<<SCB_YX_BOT_y2_array[index_z]<<" "<<SCB_YX_TOP_y2_array[index_z]<<" "<<index_z<<std::endl;

  if( eff_y/units::cm < SCB_YX_BOT_y2_array[index_z]/units::cm ) {
    double y1 = SCB_YX_BOT_y1_array[index_z];
    double x1 = SCB_YX_BOT_x1_array[index_z];
    double y2 = SCB_YX_BOT_y2_array[index_z];
    double x2 = SCB_YX_BOT_x2_array[index_z];
    double xx = (x2-x1)/(y2-y1) * (eff_y-y2) + x2;
    result = eff_x - xx;
    //std::cout<<" xp BOT "<<eff_x<<" "<<xx<<std::endl;
  }
  else if( eff_y/units::cm < SCB_YX_TOP_y2_array[index_z]/units::cm ) {
    double xx = SCB_YX_TOP_x2_array[index_z];
    result = eff_x - xx;
    //std::cout<<" xp MID "<<eff_x<<" "<<xx<<std::endl;
  }
  else {
    double y1 = SCB_YX_TOP_y1_array[index_z];
    double x1 = SCB_YX_TOP_x1_array[index_z];
    double y2 = SCB_YX_TOP_y2_array[index_z];
    double x2 = SCB_YX_TOP_x2_array[index_z];
    double xx = (x2-x1)/(y2-y1) * (eff_y-y2) + x2;
    result = eff_x - xx;
    //std::cout<<" xp TOP "<<eff_x<<" "<<xx<<std::endl;
  }

  return result;
}



bool WCPPID::ToyFiducial::inside1_outside0_SCB(WCP::Point& p, double offset_x, double tolerence_x, double tolerence_y, double tolerence_z)
{
  bool result = false;

  std::vector<double> user_SCB_YX_x_array[11], user_SCB_YX_y_array[11];
  std::vector<double> user_SCB_ZX_x_array[11], user_SCB_ZX_z_array[11];
  double space_xx = tolerence_x * units::cm;
  double space_yy = tolerence_x * units::cm;
  double space_zz = tolerence_x * units::cm;
  
  for(int idx=1; idx<=10; idx++) {
    user_SCB_YX_x_array[idx].clear();                                        user_SCB_YX_y_array[idx].clear();
    user_SCB_YX_x_array[idx].push_back(m_anode                  + space_xx); user_SCB_YX_y_array[idx].push_back(m_bottom                 + space_yy);
    user_SCB_YX_x_array[idx].push_back(SCB_YX_BOT_x1_array[idx] - space_xx); user_SCB_YX_y_array[idx].push_back(SCB_YX_BOT_y1_array[idx] + space_yy);
    user_SCB_YX_x_array[idx].push_back(SCB_YX_BOT_x2_array[idx] - space_xx); user_SCB_YX_y_array[idx].push_back(SCB_YX_BOT_y2_array[idx] + space_yy);
    user_SCB_YX_x_array[idx].push_back(SCB_YX_TOP_x2_array[idx] - space_xx); user_SCB_YX_y_array[idx].push_back(SCB_YX_TOP_y2_array[idx] - space_yy);
    user_SCB_YX_x_array[idx].push_back(SCB_YX_TOP_x1_array[idx] - space_xx); user_SCB_YX_y_array[idx].push_back(SCB_YX_TOP_y1_array[idx] - space_yy);
    user_SCB_YX_x_array[idx].push_back(m_anode                  + space_xx); user_SCB_YX_y_array[idx].push_back(m_top                    - space_yy);

    user_SCB_ZX_x_array[idx].clear();                                        user_SCB_ZX_z_array[idx].clear();
    user_SCB_ZX_x_array[idx].push_back(m_anode                  + space_xx); user_SCB_ZX_z_array[idx].push_back(m_upstream               + space_zz);
    user_SCB_ZX_x_array[idx].push_back(SCB_ZX_TOP_x1_array[idx] - space_xx); user_SCB_ZX_z_array[idx].push_back(SCB_ZX_TOP_z1_array[idx] + space_zz);
    user_SCB_ZX_x_array[idx].push_back(SCB_ZX_TOP_x2_array[idx] - space_xx); user_SCB_ZX_z_array[idx].push_back(SCB_ZX_TOP_z2_array[idx] + space_zz);
    user_SCB_ZX_x_array[idx].push_back(SCB_ZX_BOT_x2_array[idx] - space_xx); user_SCB_ZX_z_array[idx].push_back(SCB_ZX_BOT_z2_array[idx] - space_zz);
    user_SCB_ZX_x_array[idx].push_back(SCB_ZX_BOT_x1_array[idx] - space_xx); user_SCB_ZX_z_array[idx].push_back(SCB_ZX_BOT_z1_array[idx] - space_zz);
    user_SCB_ZX_x_array[idx].push_back(m_anode                  + space_xx); user_SCB_ZX_z_array[idx].push_back(m_downstream             - space_zz);   
  }
  
  
  double eff_x = p.x-offset_x;
  double eff_y = p.y;
  double eff_z = p.z;

  //std::cout<<" xp eff_x "<<eff_x<<std::endl;

  bool flag_YX = false;
  bool flag_ZX = false;

  /////////////////////////
  if( eff_z/units::cm<0 ) {
    int index_z = 1;
    flag_YX = pnpoly(user_SCB_YX_x_array[index_z], user_SCB_YX_y_array[index_z], eff_x, eff_y);
  }
  else if( eff_z/units::cm<=1037 ) {
    int index_z = (int)(eff_z/units::cm/100)+1;
    if( index_z>10 ) index_z = 10;
    flag_YX = pnpoly(user_SCB_YX_x_array[index_z], user_SCB_YX_y_array[index_z], eff_x, eff_y);
  }
  else {
    int index_z = 10;
    flag_YX = pnpoly(user_SCB_YX_x_array[index_z], user_SCB_YX_y_array[index_z], eff_x, eff_y);
  }
  
  ////////////////////////////
  if( eff_y/units::cm<-116 ) {
    int index_z = 1;
    flag_ZX = pnpoly(user_SCB_ZX_x_array[index_z], user_SCB_ZX_z_array[index_z], eff_x, eff_z);
  }
  else if( eff_y/units::cm<=117 ) {
    int index_y = (int)(eff_y/units::cm+116)/24 + 1;
    flag_ZX = pnpoly(user_SCB_ZX_x_array[index_y], user_SCB_ZX_z_array[index_y], eff_x, eff_z);
  }
  else {
    int index_z = 10;
    flag_ZX = pnpoly(user_SCB_ZX_x_array[index_z], user_SCB_ZX_z_array[index_z], eff_x, eff_z);
  }
  
  ////////////////////////////
  
  if( flag_YX && flag_ZX ) result = true;

  return result;
}


std::tuple<int, WCPPID::PR3DCluster*, WCP::Opflash*> WCPPID::ToyFiducial::M2_cosmic_tagger(double eventTime, WCP::OpflashSelection& flashes, WCPPID::PR3DCluster* main_cluster, std::vector<WCPPID::PR3DCluster*> additional_clusters, WCP::Opflash* main_flash, std::tuple<int, double, double, int>& bundle_info, WCP::Photon_Library *pl, int time_offset, int nrebin, float unit_dis, WCP::ToyCTPointCloud& ct_point_cloud, int run_no, int subrun_no, int event_no, bool flag_data, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map, bool flag_timestamp, bool debug_tagger)
{
  std::cout<<std::endl<<" ---------> Hello M2"<<std::endl<<std::endl;
  //std::cout<<" M2: main_flash->get_time() "<<main_flash->get_time()<<std::endl;

  int flag_M2 = 0;// 0 pass, 1 tgm with flash, 2 tgm wo flash, 3 stm candidates

  WCP::Opflash* new_flash = main_flash;
    
  /// debug  

  double meas_posZ_center = 0;
  int num_meas_posZ = 0;
  double pmt_posZ_center = 0;

  {
    ToyPointCloud *pcloud = main_cluster->get_point_cloud();
    if (pcloud!=0){
      WCP::WCPointCloud<double>& cloud = pcloud->get_cloud();
      for (size_t i=0;i!=cloud.pts.size();i++){
	double x = cloud.pts[i].x/units::cm;
	double y = cloud.pts[i].y/units::cm;
	double z = cloud.pts[i].z/units::cm;

	num_meas_posZ++;
	meas_posZ_center += z;
      }
    }
    meas_posZ_center = meas_posZ_center/num_meas_posZ;

    std::cout<<TString::Format(" ccheck %6d %6d %6d, Zcenter %8.4f", run_no, subrun_no, event_no, meas_posZ_center)<<std::endl;
  }
  ///

  double pmt_posZ_cm_left2right[21] = {0,
				       51.37,  85.61, 125.57, 172.66, 
				       241.15, 283.96, 325.35, 371.01,
				       448.07, 496.58, 535.11, 579.35,
				       660.69, 704.92, 746.31, 789.12,
				       859.04, 904.70, 941.80, 981.76};
  
  // pmt map, DocDB-7555
  double pmt_posZ_cm_array[33] = {0};

  pmt_posZ_cm_array[1] = pmt_posZ_cm_left2right[2];
  pmt_posZ_cm_array[2] = pmt_posZ_cm_left2right[3];
  pmt_posZ_cm_array[3] = pmt_posZ_cm_left2right[1];
  pmt_posZ_cm_array[4] = pmt_posZ_cm_left2right[4];
  pmt_posZ_cm_array[5] = pmt_posZ_cm_left2right[1];
  pmt_posZ_cm_array[6] = pmt_posZ_cm_left2right[2];
  pmt_posZ_cm_array[7] = pmt_posZ_cm_left2right[3];
  
  pmt_posZ_cm_array[8]  = pmt_posZ_cm_left2right[6];
  pmt_posZ_cm_array[9]  = pmt_posZ_cm_left2right[7];
  pmt_posZ_cm_array[10] = pmt_posZ_cm_left2right[5];
  pmt_posZ_cm_array[11] = pmt_posZ_cm_left2right[8];
  pmt_posZ_cm_array[12] = pmt_posZ_cm_left2right[6];
  pmt_posZ_cm_array[13] = pmt_posZ_cm_left2right[7];

  pmt_posZ_cm_array[14] = pmt_posZ_cm_left2right[10];
  pmt_posZ_cm_array[15] = pmt_posZ_cm_left2right[11];
  pmt_posZ_cm_array[16] = pmt_posZ_cm_left2right[9];
  pmt_posZ_cm_array[17] = pmt_posZ_cm_left2right[12];
  pmt_posZ_cm_array[18] = pmt_posZ_cm_left2right[10];
  pmt_posZ_cm_array[19] = pmt_posZ_cm_left2right[11];

  pmt_posZ_cm_array[20] = pmt_posZ_cm_left2right[14];
  pmt_posZ_cm_array[21] = pmt_posZ_cm_left2right[15];
  pmt_posZ_cm_array[22] = pmt_posZ_cm_left2right[13];
  pmt_posZ_cm_array[23] = pmt_posZ_cm_left2right[16];
  pmt_posZ_cm_array[24] = pmt_posZ_cm_left2right[14];
  pmt_posZ_cm_array[25] = pmt_posZ_cm_left2right[15];

  pmt_posZ_cm_array[26] = pmt_posZ_cm_left2right[18];
  pmt_posZ_cm_array[27] = pmt_posZ_cm_left2right[19];
  pmt_posZ_cm_array[28] = pmt_posZ_cm_left2right[20];
  pmt_posZ_cm_array[29] = pmt_posZ_cm_left2right[17];
  pmt_posZ_cm_array[30] = pmt_posZ_cm_left2right[20];
  pmt_posZ_cm_array[31] = pmt_posZ_cm_left2right[18];
  pmt_posZ_cm_array[32] = pmt_posZ_cm_left2right[19];
  


  //Take in flash and cluster info
  std::vector<std::vector<WCPointCloud<double>::WCPoint>> extreme_points = main_cluster->get_extreme_wcps();
  double time_slice_width = nrebin * unit_dis * 0.5 * units::mm;
  double main_offset_x = (main_flash->get_time() - time_offset)*2./nrebin*time_slice_width;

  //Check whether the track is at least 10 cm long
  Point p0(extreme_points[0][0].x,extreme_points[0][0].y,extreme_points[0][0].z);
  Point p1(extreme_points[1][0].x,extreme_points[1][0].y,extreme_points[1][0].z);

  //std::cout<<" xpdebug "<<p0.x/units::cm<<" "<<p0.y/units::cm<<" "<<p0.z/units::cm<<" "<<1<<std::endl;
  //std::cout<<" xpdebug "<<p1.x/units::cm<<" "<<p1.y/units::cm<<" "<<p1.z/units::cm<<" "<<1<<std::endl;


  Double_t cos_pe_low[32]={11,11,11,11,10,
                           7,8,8,10,7,
                           11,11,11,11,10,
                           9,11,11,7,9,
                           11,10,11,11,11,
                           11,11,10,11,11,
                           9,10};
  Double_t cos_pe_mid[32]={34,32,28,35,22,
                           23,22,24,33,30,
                           35,35,33,36,33,
                           33,36,33,19,27,
                           32,23,42,32,33,
                           34,34,24,33,35,
                           25,32};


  double distance = sqrt(pow(p0.x-p1.x,2)+pow(p0.y-p1.y,2)+pow(p0.z-p1.z,2));
  if(distance > 10*units::cm){
    //iterate over all the flashes to check each flash for STM/TGM conditions
    
    double bundle_ks = std::get<1>(bundle_info);
    double bundle_chi2 = std::get<2>(bundle_info);
    int bundle_ndf = std::get<3>(bundle_info);
    std::cout<<"xpks "<<bundle_ks<<"\t"<<bundle_chi2<<"\t"<<bundle_ndf<<std::endl;

    double diff_pe_meas_pred_origin = 0;
    double origin_meas_pe_total = 0;
    double origin_pred_pe_total = 0;
    for (auto it = flashes.begin(); it!= flashes.end(); it++){
      Opflash *flash = (*it);
      double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
      double flash_time = flash->get_time();

      if(flash->get_time() == main_flash->get_time()){	

	std::vector<double> pred_pmt_light = calculate_pred_pe(eventTime, run_no, offset_x, pl, main_cluster, additional_clusters, flash, flag_data, flag_timestamp);	
	for(int i=0; i<int(pred_pmt_light.size()); i++){
	  double flash_pe = flash->get_PE(i);
	  double pred_pe = pred_pmt_light[i];
	  diff_pe_meas_pred_origin += fabs( flash_pe - pred_pe );
	  origin_meas_pe_total += flash_pe;
	  origin_pred_pe_total += pred_pe;	  
	  //std::cout<<" origin ---> "<<i<<"\t"<<flash_pe<<"\t"<<pred_pmt_light[i]<<"\t"<<diff_pe_meas_pred_origin<<std::endl;

	  //std::cout<<TString::Format(" ddcheck %6d %6d %6d pmt %2d pe %8.2f",
	  //			     run_no, subrun_no, event_no, i+1, flash_pe
	  //			     )<<std::endl;

	  pmt_posZ_center += pmt_posZ_cm_array[i+1] * flash_pe;
	}
	
	pmt_posZ_center = pmt_posZ_center/origin_meas_pe_total;
	pmt_posZ_center = 1037 - pmt_posZ_center;

	std::cout<<TString::Format(" ddcheck %6d %6d %6d, Zcenter %8.4f %8.4f", 
				   run_no, subrun_no, event_no, meas_posZ_center, pmt_posZ_center)<<std::endl;
	
	////////////////////////////////////////////////////////////////////////////////

	TH1D *h1_meas = new TH1D("h1_meas", "h1_meas", 32, 0, 32);
	TH1D *h1_pred = new TH1D("h1_pred", "h1_pred", 32, 0, 32);
    
	double chi2_user = 0;
	int ndf_user = 0;

	for(int i=0; i<int(pred_pmt_light.size()); i++){
	  double flash_pe = flash->get_PE(i);
	  double pred_pe = pred_pmt_light[i];
	  double pe_err = flash->get_PE_err(i);

	  int j = i;
	  if( (pred_pe < cos_pe_low[j] || (pred_pe < cos_pe_mid[j]*1.1 && flash_pe==0 ) ) && flash->get_type()==1 ) {
	    pred_pe = 0;
	  }	    
	    
	  h1_meas->SetBinContent(i+1, flash_pe);
	  h1_pred->SetBinContent(i+1, pred_pe);
	  //std::cout<<" ---> "<<i+1<<"\t"<<flash_pe<<"\t"<<pred_pe<<std::endl;

	  ///
	  chi2_user += pow(pred_pe - flash_pe,2)/pow(pe_err,2);
	  if( pred_pe==0 && flash_pe==0 ) {
	  }
	  else {
	    ndf_user++;
	  }

	}

	double temp_ks_dis = 1;
	if( h1_pred->GetSum()!=0 ) {
	  temp_ks_dis = h1_meas->KolmogorovTest(h1_pred,"M");
	}
	  
	std::cout<<" checkxp original ks "<<temp_ks_dis<<"\t"<<"chi "<<chi2_user<<"\t"<<ndf_user<<"\t ndf "<<bundle_ndf<<"\t"<<ndf_user<<std::endl;
	  
	delete h1_meas;
	delete h1_pred;

	break;
      }// if(flash->get_time() == main_flash->get_time())

    }
    

    ///////////////// tgm with flash
    /////////////////

    TH1D *h1_meas = new TH1D("h1_meas", "h1_meas", 32, 0, 32);
    TH1D *h1_pred = new TH1D("h1_pred", "h1_pred", 32, 0, 32);
    
    for (auto it = flashes.begin(); it!= flashes.end(); it++){
      Opflash *flash = (*it);
      double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
      double flash_time = flash->get_time();

      //Skip the already matched flash
      if(flash->get_time() == main_flash->get_time()){
	continue;
      }

      double PE_ratio_cut = 0.8;
      double dist2anode_cut = 2*units::cm;
      double original_PE_KS_cut = 0.05;// GE 0.1 --> not well matched between meas and pred

      if( bundle_ks < original_PE_KS_cut ) {// check original ks: current ks is good
	break;
      }

      h1_meas->Reset();
      h1_pred->Reset();

      double tolerence_x = 2;// cm
      double tolerence_y = 2;// cm
      double tolerence_z = 2;// cm

      bool flag_small_SCB_p0 = false;
      bool flag_large_SCB_p0 = false;
      bool flag_small_SCB_p1 = false;
      bool flag_large_SCB_p1 = false;

      flag_small_SCB_p0 = inside1_outside0_SCB(p0, offset_x, tolerence_x, tolerence_y, tolerence_z);
      flag_large_SCB_p0 = inside1_outside0_SCB(p0, offset_x, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1));
      flag_small_SCB_p1 = inside1_outside0_SCB(p1, offset_x, tolerence_x, tolerence_y, tolerence_z);
      flag_large_SCB_p1 = inside1_outside0_SCB(p1, offset_x, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1));
      

      if( flag_small_SCB_p0==false && flag_small_SCB_p1==false && flag_large_SCB_p0==true && flag_large_SCB_p1==true ) {// check boundary
		
	////
	bool flag_NotNear_anode = true;
	if( (p0.x-offset_x)<dist2anode_cut || (p1.x-offset_x)<dist2anode_cut ) flag_NotNear_anode = false;
	
	if( flag_NotNear_anode ) {

	  double chi2_user = 0;
	  int ndf_user = 0;

	  std::vector<double> pred_pmt_light = calculate_pred_pe(eventTime, run_no, offset_x, pl, main_cluster, additional_clusters, flash, flag_data, flag_timestamp);	
	  for(int i=0; i<int(pred_pmt_light.size()); i++){
	    double flash_pe = flash->get_PE(i);
	    double pred_pe = pred_pmt_light[i];
	    double pe_err = flash->get_PE_err(i);

	    int j = i;
	    if( (pred_pe < cos_pe_low[j] || (pred_pe < cos_pe_mid[j]*1.1 && flash_pe==0 ) ) && flash->get_type()==1 ) {
	      pred_pe = 0;
	    }	    
	    
	    h1_meas->SetBinContent(i+1, flash_pe);
	    h1_pred->SetBinContent(i+1, pred_pe);
	    //std::cout<<" ---> "<<i+1<<"\t"<<flash_pe<<"\t"<<pred_pe<<std::endl;

	    ///
	    chi2_user += pow(pred_pe - flash_pe,2)/pow(pe_err,2);
	    if( pred_pe==0 && flash_pe==0 ) {
	    }
	    else {
	      ndf_user++;
	    }
	    
	  }

	  double temp_ks_dis = 1;
	  if( h1_pred->GetSum()!=0 ) {
	    temp_ks_dis = h1_meas->KolmogorovTest(h1_pred,"M");
	  }

	  std::cout<<TString::Format(" checkxp pass: bounary + original_ks + not near anode, ks org/new %6.4f %6.4f, chi2 org/new %8.2f %8.2f, ndf %2d %2d, %7d %7d %7d", 
				     bundle_ks, temp_ks_dis,
				     bundle_chi2, chi2_user,
				     bundle_ndf, ndf_user, 
				     run_no, subrun_no, event_no
				     )<<std::endl;
	  
	  if( chi2_user < bundle_chi2 ) {

	    
	    
	    
	    bool flag_check_tgm = check_tgm(main_cluster,flash,offset_x,ct_point_cloud);
	    
	    if( flag_check_tgm ) {
	      new_flash = flash;
	      
	      int flag_M2_middle = 1;
	      flag_M2 = 1;
	      std::cout<<" M2 "<<flag_M2_middle<<"\t"<<flag_M2<<"\t"<<run_no<<"\t"<<subrun_no<<"\t"<<event_no
		       <<" tgm wi flash, "<<std::endl;
	      
	      std::cout<<TString::Format(" xxcheck original chi2/ndf %8.3f %2d, new %8.3f %2d, run %6d %6d %6d", 
					 bundle_chi2, bundle_ndf, 
					 chi2_user, ndf_user,
					 run_no, subrun_no, event_no
					 )<<std::endl;

	      break;
	    }// tgm_tagger
	    
	  }// chi2_current < chi2_original

	}// Near anode
	
      }// check boundary

    }// loop flash

    delete h1_meas;
    delete h1_pred;

    ///////////////// tgm without flash
    /////////////////

    if( flag_M2==0 ) {
      double tolerence_x = 2;// cm
      double tolerence_y = 2;// cm
      double tolerence_z = 2;// cm
      
      double original_PE_KS_cut = 0.1;// GE 0.1 --> not well matched between meas and pred

      if( bundle_ks>=original_PE_KS_cut ) {// original ks is not good

	if( flag_M2==0 ) {// case1
	  bool flag_small_SCB_p0 = false;
	  bool flag_large_SCB_p0 = false;
	  bool flag_small_SCB_p1 = false;
	  bool flag_large_SCB_p1 = false;
	
	  double user_offset = M2_offset_YX_x( p0 );
		
	  flag_small_SCB_p0 = inside1_outside0_SCB(p0, user_offset, tolerence_x, tolerence_y, tolerence_z);
	  flag_large_SCB_p0 = inside1_outside0_SCB(p0, user_offset, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1));
	  flag_small_SCB_p1 = inside1_outside0_SCB(p1, user_offset, tolerence_x, tolerence_y, tolerence_z);
	  flag_large_SCB_p1 = inside1_outside0_SCB(p1, user_offset, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1));
	
	  if( flag_small_SCB_p0==false && flag_small_SCB_p1==false && flag_large_SCB_p0==true && flag_large_SCB_p1==true ) {// check boundary

	    bool flag_tgm_tagger = M2_check_tgm( main_cluster, user_offset, ct_point_cloud );

	    if( flag_tgm_tagger ) {

	      bool flag_AllExtreeamP_inside = true;
	    
	      ToyPointCloud *pcloud = main_cluster->get_point_cloud();
	      if (pcloud!=0){
		WCP::WCPointCloud<double>& cloud = pcloud->get_cloud();
		for (size_t i=0;i!=cloud.pts.size();i++){
		  double x = cloud.pts[i].x;
		  double y = cloud.pts[i].y;
		  double z = cloud.pts[i].z;
		  Point pobj(x,y,z);
		  if( inside1_outside0_SCB(pobj, user_offset, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1))==false ) {
		    flag_AllExtreeamP_inside = false;
		    break;
		  }
		}
	      }
	     
	      if( flag_AllExtreeamP_inside ) {

		int flag_M2_middle = 2;
		flag_M2 = 2;
		new_flash = NULL;

		std::cout<<" M2 "<<flag_M2_middle<<"\t"<<flag_M2<<"\t"<<run_no<<"\t"<<subrun_no<<"\t"<<event_no<<"\t"
			 <<" case1, "
			 <<std::endl;

	      }// flag_AllExtreeamP_inside
	    }// tgm tagger
	  }// check boundary
	}// case1


	if( flag_M2==0 ) {// case2
	  bool flag_small_SCB_p0 = false;
	  bool flag_large_SCB_p0 = false;
	  bool flag_small_SCB_p1 = false;
	  bool flag_large_SCB_p1 = false;
	
	  double user_offset = M2_offset_YX_x( p1 );
		
	  flag_small_SCB_p0 = inside1_outside0_SCB(p0, user_offset, tolerence_x, tolerence_y, tolerence_z);
	  flag_large_SCB_p0 = inside1_outside0_SCB(p0, user_offset, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1));
	  flag_small_SCB_p1 = inside1_outside0_SCB(p1, user_offset, tolerence_x, tolerence_y, tolerence_z);
	  flag_large_SCB_p1 = inside1_outside0_SCB(p1, user_offset, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1));
	
	  if( flag_small_SCB_p0==false && flag_small_SCB_p1==false && flag_large_SCB_p0==true && flag_large_SCB_p1==true ) {// check boundary

	    bool flag_tgm_tagger = M2_check_tgm( main_cluster, user_offset, ct_point_cloud );

	    if( flag_tgm_tagger ) {

	      bool flag_AllExtreeamP_inside = true;
	    
	      ToyPointCloud *pcloud = main_cluster->get_point_cloud();
	      if (pcloud!=0){
		WCP::WCPointCloud<double>& cloud = pcloud->get_cloud();
		for (size_t i=0;i!=cloud.pts.size();i++){
		  double x = cloud.pts[i].x;
		  double y = cloud.pts[i].y;
		  double z = cloud.pts[i].z;
		  Point pobj(x,y,z);
		  if( inside1_outside0_SCB(pobj, user_offset, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1))==false ) {
		    flag_AllExtreeamP_inside = false;
		    break;
		  }
		}
	      }
	     
	      if( flag_AllExtreeamP_inside ) {

		int flag_M2_middle = 2;
		flag_M2 = 2;
		new_flash = NULL;

		std::cout<<" M2 "<<flag_M2_middle<<"\t"<<flag_M2<<"\t"<<run_no<<"\t"<<subrun_no<<"\t"<<event_no<<"\t"
			 <<" case2, "
			 <<std::endl;

	      }// flag_AllExtreeamP_inside
	    }// tgm tagger
	  }// check boundary
	}// case2

      }// original ks is not good

    }

    ///////////////// stm candidates
    /////////////////

    if( flag_M2==0 ) {
      
      double tolerence_x = 2;// cm
      double tolerence_y = 2;// cm
      double tolerence_z = 2;// cm
      	  
      double dist2anode_cut = 2*units::cm;
      double original_PE_KS_cut = 0.05;// GE 0.1 --> not well matched between meas and pred

      for (auto it = flashes.begin(); it!= flashes.end(); it++){
	Opflash *flash = (*it);
	double offset_x = (flash->get_time() - time_offset)*2./nrebin*time_slice_width;
	double flash_time = flash->get_time();

	//Skip the already matched flash
	if(flash->get_time() == main_flash->get_time()){
	  continue;
	}

	if( bundle_ks < original_PE_KS_cut ) {// check original ks: original ks is good
	  break;
	}

	///////////////////

	bool flag_NotNear_anode = true;
	if( (p0.x-offset_x)<dist2anode_cut || (p1.x-offset_x)<dist2anode_cut ) flag_NotNear_anode = false;

	if( flag_NotNear_anode ) {
	
	  double chi2_user = 0;
	  int ndf_user = 0;

	  double pmt_posZ_center = 0;
	  double total_meas_pe = 0;

	  std::vector<double> pred_pmt_light = calculate_pred_pe(eventTime, run_no, offset_x, pl, main_cluster, additional_clusters, flash, flag_data, flag_timestamp);	
	  for(int i=0; i<int(pred_pmt_light.size()); i++){
	    double flash_pe = flash->get_PE(i);
	    double pred_pe = pred_pmt_light[i];
	    double pe_err = flash->get_PE_err(i);

	    int j = i;
	    if( (pred_pe < cos_pe_low[j] || (pred_pe < cos_pe_mid[j]*1.1 && flash_pe==0 ) ) && flash->get_type()==1 ) {
	      pred_pe = 0;
	    }	    

	    ///
	    chi2_user += pow(pred_pe - flash_pe,2)/pow(pe_err,2);
	    if( pred_pe==0 && flash_pe==0 ) {
	    }
	    else {
	      ndf_user++;
	    }
	    
	    pmt_posZ_center += pmt_posZ_cm_array[i+1] * flash_pe;
	    total_meas_pe += flash_pe;
	    
	  }// loop pmt

	  pmt_posZ_center = pmt_posZ_center/total_meas_pe;
	  pmt_posZ_center = 1037 - pmt_posZ_center;

	  if( chi2_user < bundle_chi2 && fabs(pmt_posZ_center-meas_posZ_center)<100 ) {
	  
	    for(int i=0; i<=1; i++) {
	      for(int j=0; j<=0; j++) {

		bool flag_small_SCB_p0 = false;
		bool flag_large_SCB_p0 = false;
		bool flag_small_SCB_p1 = false;
		bool flag_large_SCB_p1 = false;
	    
		Point pobj(extreme_points[i][j].x,extreme_points[i][j].y,extreme_points[i][j].z);	  
	    
		flag_small_SCB_p0 = inside1_outside0_SCB(pobj, offset_x, tolerence_x, tolerence_y, tolerence_z);
		flag_large_SCB_p0 = inside1_outside0_SCB(pobj, offset_x, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1));
	    
		if( flag_small_SCB_p0==false && flag_large_SCB_p0==true ) {
		  
		  bool flag_AllExtreeamP_inside = true;
		  
		  ToyPointCloud *pcloud = main_cluster->get_point_cloud();
		  if (pcloud!=0){
		    WCP::WCPointCloud<double>& cloud = pcloud->get_cloud();
		    for (size_t i=0;i!=cloud.pts.size();i++){
		      double x = cloud.pts[i].x;
		      double y = cloud.pts[i].y;
		      double z = cloud.pts[i].z;
		      Point pobj(x,y,z);
		      if( inside1_outside0_SCB(pobj, offset_x, tolerence_x*(-1), tolerence_y*(-1), tolerence_z*(-1))==false ) {
			flag_AllExtreeamP_inside = false;
			break;
		      }
		    }
		  }

		  if( flag_AllExtreeamP_inside ) {
		    bool flag_user_stm =  M2_check_stm(main_cluster, additional_clusters, offset_x, flash_time, ct_point_cloud, global_wc_map);
		  
		    if( flag_user_stm ) {
		      int flag_M2_middle = 3;
		      flag_M2 = 3;
		      new_flash = flash;
		  
		      std::cout<<" M2 "<<flag_M2_middle<<"\t"<<flag_M2<<"\t"<<run_no<<"\t"<<subrun_no<<"\t"<<event_no <<"\t stmxp3 "
			       <<TString::Format(" ks %6.4f", bundle_ks )
			       <<std::endl;		  
		      if( flag_M2!=0 ) break;
		    }// stm tagger
		  }// flag_AllExtreeamP_inside
		}// check boundary
	      }// loop j
	      if( flag_M2!=0 ) break;
	    }// loop i
	  }// new chi2 < original chi2
	}// flag_NotNear_anode
	if( flag_M2!=0 ) break;
      }// loop flash
    }// if( flag_M2==0 )



  }// if(distance > 10*units::cm)


  int tag_type = flag_M2;// 0-normal, 1-tgm wi flash, 2-tgm wo flash, 3-stm wi flash

  return std::make_tuple(tag_type, main_cluster, new_flash);
}


bool WCPPID::ToyFiducial::M2_check_tgm(WCPPID::PR3DCluster* main_cluster, double offset_x, WCP::ToyCTPointCloud& ct_point_cloud) {

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
 
      if ((!flag_p1_inside) && (!flag_p2_inside)){

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
	  if (1){
		    
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
	    if (1){
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

	}
      }
    }
  }
  

  return false;

}





bool WCPPID::ToyFiducial::M2_check_stm(WCPPID::PR3DCluster* main_cluster, std::vector<WCPPID::PR3DCluster*>& additional_clusters, double offset_x, double flash_time, WCP::ToyCTPointCloud& ct_point_cloud, std::map<int,std::map<const WCP::GeomWire*, WCP::SMGCSelection > >& global_wc_map){

  //  check_full_detector_dead();
  //std::cout << flag_other_clusters << std::endl;
  
  TVector3 drift_dir(1,0,0);
  // hard coded for U and V plane ... 
  TVector3 U_dir(0,cos(60./180.*3.1415926),sin(60./180.*3.1415926));
  TVector3 V_dir(0,cos(60./180.*3.1415926),-sin(60./180.*3.1415926));
  TVector3 W_dir(0,1,0);

  Vector main_dir = main_cluster->get_PCA_axis(0);
  TVector3 dir_main(main_dir.x,main_dir.y,main_dir.z);
  
  std::vector<WCPointCloud<double>::WCPoint> candidate_exit_wcps;
  std::set<int> temp_set;
  std::pair<WCPointCloud<double>::WCPoint,WCPointCloud<double>::WCPoint> wcps;
  
  // first round check
  if (main_cluster->get_point_cloud_steiner()->get_cloud().pts.size()==0)
    return false;

  {
    // get two extreme points ...
    wcps = main_cluster->get_two_boundary_wcps(2,true);
    // figure out the end points ...
    std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster->get_extreme_wcps();
    
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.first);
      out_vec_wcps.push_back(temp_wcps);
    }
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.second);
      out_vec_wcps.push_back(temp_wcps);
    }
    
    // boundary check
    for (size_t i=0;i!=out_vec_wcps.size();i++){
      bool flag_save = false;
      // check all the points ... 
      for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){

	Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
	//std::cout << p1 << " " << inside_fiducial_volume(p1,offset_x) << std::endl;
	if (!inside_fiducial_volume(p1,offset_x)){
	  candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  flag_save = true;
	  break;
	}
      }
      
      if (!flag_save){
	// check direction
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
	
	
	if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	  if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
	    flag_save = true;
	    candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  }
	}
	
	if (!flag_save){
	  if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
	    if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	      flag_save = true;
	      candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	    }
	  }
	}
      }
    }
    
    //std::cout << wcps.first.x << " " <<wcps.first.y << " " << wcps.first.z << std::endl;
    //std::cout << wcps.second.x << " " <<wcps.second.y << " " << wcps.second.z << std::endl;
    
    for (size_t i=0; i!=candidate_exit_wcps.size(); i++){
      double dis1 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.first.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.first.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.first.z,2));
      double dis2 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.second.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.second.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.second.z,2));
      
      //      std::cout << candidate_exit_wcps.at(i).x << " " << candidate_exit_wcps.at(i).y << " " << candidate_exit_wcps.at(i).z << " " << dis1 << " " << dis2 << std::endl;
      
      // essentially one of the extreme points ...
      if (dis1 < dis2){
	if (dis1 < 1.0*units::cm)  temp_set.insert(0);
      }else{
	if (dis2 < 1.0*units::cm)  temp_set.insert(1);
      }
    }

    //    std::cout << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;
    
    // protection against two end point situation
    if (temp_set.size()==2){
      WCP::Point tp1(wcps.first.x,wcps.first.y,wcps.first.z);
      WCP::Point tp2(wcps.second.x,wcps.second.y,wcps.second.z);
      
      temp_set.clear();
      
      if ((!inside_fiducial_volume(tp1,offset_x))) temp_set.insert(0);
      if ((!inside_fiducial_volume(tp2,offset_x))) temp_set.insert(1);
      if (temp_set.size()==0){
	temp_set.insert(0);
	temp_set.insert(1);
      }
    }
  }

  if (temp_set.size()==0){
    candidate_exit_wcps.clear();
    // get two extreme points ...
    wcps = main_cluster->get_two_boundary_wcps(2);
    // figure out the end points ...
    std::vector<std::vector<WCPointCloud<double>::WCPoint>> out_vec_wcps = main_cluster->get_extreme_wcps();
    
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.first);
      out_vec_wcps.push_back(temp_wcps);
    }
    {
      std::vector<WCPointCloud<double>::WCPoint> temp_wcps;
      temp_wcps.push_back(wcps.second);
      out_vec_wcps.push_back(temp_wcps);
    }
    
    // boundary check
    for (size_t i=0;i!=out_vec_wcps.size();i++){
      bool flag_save = false;
      // check all the points ... 
      for (size_t j=0;j!=out_vec_wcps.at(i).size();j++){
	Point p1(out_vec_wcps.at(i).at(j).x,out_vec_wcps.at(i).at(j).y,out_vec_wcps.at(i).at(j).z);
	if (!inside_fiducial_volume(p1,offset_x)){
	  candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  flag_save = true;
	  break;
	}
      }
      
      if (!flag_save){
	// check direction
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
	
	
	if ( (angle1_1 < 10 || angle2_1 < 10 || angle3_1 < 5)){
	  if (!check_signal_processing(p1,dir,ct_point_cloud,1*units::cm,offset_x)){
	    flag_save = true;
	    candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	  }
	}

	//	std::cout << fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.) << " " << check_dead_volume(p1,dir,1*units::cm,offset_x) << std::endl;
	if (!flag_save){
	  if (fabs((3.1415926/2.-dir.Angle(dir_main))/3.1415926*180.)>60 ){
	    if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	      flag_save = true;
	      candidate_exit_wcps.push_back(out_vec_wcps.at(i).at(0));
	    }
	  }
	}
      }
    }
    
    //std::cout << wcps.first.x << " " <<wcps.first.y << " " << wcps.first.z << std::endl;
    //std::cout << wcps.second.x << " " <<wcps.second.y << " " << wcps.second.z << std::endl;
    
    for (size_t i=0; i!=candidate_exit_wcps.size(); i++){
      double dis1 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.first.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.first.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.first.z,2));
      double dis2 = sqrt(pow(candidate_exit_wcps.at(i).x - wcps.second.x,2) + pow(candidate_exit_wcps.at(i).y - wcps.second.y,2) + pow(candidate_exit_wcps.at(i).z - wcps.second.z,2));
      
      //std::cout << candidate_exit_wcps.at(i).x << " " << candidate_exit_wcps.at(i).y << " " << candidate_exit_wcps.at(i).z << " " << dis1 << " " << dis2 << std::endl;
      
      // essentially one of the extreme points ...
      if (dis1 < dis2){
	if (dis1 < 1.0*units::cm)  temp_set.insert(0);
      }else{
	if (dis2 < 1.0*units::cm)  temp_set.insert(1);
      }
    }

    //std::cout << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;
    
    
    // protection against two end point situation
    if (temp_set.size()==2){
      WCP::Point tp1(wcps.first.x,wcps.first.y,wcps.first.z);
      WCP::Point tp2(wcps.second.x,wcps.second.y,wcps.second.z);
      
      temp_set.clear();
      
      if ((!inside_fiducial_volume(tp1,offset_x))) temp_set.insert(0);
      if ((!inside_fiducial_volume(tp2,offset_x))) temp_set.insert(1);
      if (temp_set.size()==0){
	temp_set.insert(0);
	temp_set.insert(1);
      }
    }
  }
  

  
  // fully contained, so not a STM
  if (candidate_exit_wcps.size()==0) {
    std::cout << "Mid Point: A" << std::endl;
    return false;
  }
  
  std::cout << "end_point: " << temp_set.size() << " " << candidate_exit_wcps.size() << std::endl;

  // It is possible that we have two points out of fiducial
  // Michel electron is outside the boundary ...

  
  
  // Crawl backward according to the graph ??? ...
  WCPointCloud<double>::WCPoint first_wcp;
  WCPointCloud<double>::WCPoint last_wcp;
  bool flag_double_end = false;

  if (temp_set.size()!=0){
    if (*temp_set.begin()==0){
      first_wcp = wcps.first;
      last_wcp = wcps.second;
    }else{
      first_wcp = wcps.second;
      last_wcp = wcps.first;
    }
    if (temp_set.size()==2) flag_double_end = true;
    
  }else{
    if (candidate_exit_wcps.size()==1){
      first_wcp = candidate_exit_wcps.at(0);
      TVector3 dir1(wcps.first.x - candidate_exit_wcps.at(0).x,
		    wcps.first.y - candidate_exit_wcps.at(0).y,
		    wcps.first.z - candidate_exit_wcps.at(0).z);
      TVector3 dir2(wcps.second.x - candidate_exit_wcps.at(0).x,
		    wcps.second.y - candidate_exit_wcps.at(0).y,
		    wcps.second.z - candidate_exit_wcps.at(0).z);
      double dis1 = dir1.Mag();
      double dis2 = dir2.Mag();

      if (dir1.Angle(dir2) > 120/180.*3.1415926 && dis1 > 20*units::cm &&
      	  dis2 > 20*units::cm){
	std::cout << "Mid Point: B" << std::endl;
	return false;
      }else{
	if (dis1 < dis2){
	  last_wcp = wcps.second;	   
	}else{
	  last_wcp = wcps.first; 
	}
      }
      
    }else{
      std::cout << "Mid Point: C" << std::endl;
      return false;
    }
  }

 
  bool flag_other_clusters = check_other_clusters(main_cluster, additional_clusters);
  //  std::cout << "haha " << flag_other_clusters << std::endl;
  
  // forward check ...
  {
    if (flag_double_end) std::cout << "Forward check! " << std::endl;
    // regular crawling ...
    main_cluster->do_rough_path(first_wcp, last_wcp);
    // std::cout << "haha" << std::endl;
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    //std::cout << "haha" << std::endl;
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
    //std::cout << "test" << std::endl; 
    if (main_cluster->get_fine_tracking_path().size()<=3) return false;
    
    //std::cout << "haha " << std::endl;
    Point mid_p = main_cluster->adjust_rough_path(); 
    //std::cout << "haha" << std::endl;
    // fitting trajectory and dQ/dx...
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);

    //    std::cout << "test" << std::endl; 
    // check both end points for TGM ...
    WCP::PointVector& pts = main_cluster->get_fine_tracking_path();
    std::vector<double>& dQ = main_cluster->get_dQ();
    std::vector<double>& dx = main_cluster->get_dx();
    
    int kink_num = find_first_kink(main_cluster);
    
    double left_L = 0; 
    double left_Q = 0;
    double exit_L = 0; 
    double exit_Q = 0;
    for (size_t i=0;i!=kink_num;i++){
      exit_L += dx.at(i);
      exit_Q += dQ.at(i);
    }
    for (size_t i = kink_num; i!=dx.size(); i++){
      left_L += dx.at(i);
      left_Q += dQ.at(i);
    }
    
    std::cout << "Left: " << exit_L/units::cm << " " << left_L/units::cm << " " << (left_Q/(left_L/units::cm+1e-9))/50e3 << " " << (exit_Q/(exit_L/units::cm+1e-9)/50e3) << std::endl;

    //std::cout << pts.front() << " " << pts.back() << " " << (!inside_fiducial_volume(pts.front(),offset_x)) << " " << (!inside_fiducial_volume(pts.back(),offset_x)) << std::endl;
    if ( (!inside_fiducial_volume(pts.front(),offset_x)) && (!inside_fiducial_volume(pts.back(),offset_x))){
      // std::cout << exit_L_TGM/units::cm << " " << left_L_TGM/units::cm << std::endl;
      bool flag_TGM_anode = false;
      
      if ((pts.back().x < 2*units::cm || pts.front().x < 2*units::cm) && kink_num >=0 && kink_num < pts.size()){ // at Anode ...
	if (pts.at(kink_num).x < 6*units::cm){
	  TVector3 v10(pts.back().x-pts.at(kink_num).x,pts.back().y-pts.at(kink_num).y,pts.back().z-pts.at(kink_num).z);
	  TVector3 v20(pts.front().x-pts.at(kink_num).x,pts.front().y-pts.at(kink_num).y,pts.front().z-pts.at(kink_num).z);
	  if (fabs(v10.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v10.Mag()>15*units::cm || fabs(v20.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v20.Mag()>15*units::cm)
	    flag_TGM_anode = true;
	  //std::cout << v10.Angle(drift_dir)/3.1415926*180. << " " << v20.Angle(drift_dir)/3.1415926*180. << std::endl;
	}
      }
      if ((exit_L < 3*units::cm || left_L < 3*units::cm || flag_TGM_anode)){
	std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
	//event_type |= 1UL << 3;
	return false;
      }
    }else if ((!inside_fiducial_volume(pts.front(),offset_x)) && left_L < 3*units::cm){
      // check dead volume and SP ...
      Point p1 = pts.back();
      TVector3 dir = main_cluster->VHoughTrans(p1,30*units::cm);
      dir *= (-1);
      if (!check_dead_volume(p1,dir,1*units::cm,offset_x)){
	if (exit_L < 3*units::cm || left_L < 3*units::cm){
	  std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
	  //event_type |= 1UL << 3;
	  return false;
	}
      }
      //std::cout << "ABC: " << dir.X() << " " << dir.Y() << " " << dir.Z() << " " << check_dead_volume(p1,dir,1*units::cm,offset_x) << std::endl;
    }
  
    if (left_L > 40*units::cm || left_L > 7.5*units::cm && (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0){
      if (!flag_double_end){
	std::cout << "Mid Point A " << inside_dead_region(mid_p) << " " << mid_p << " " << left_L  << " " << (left_Q/(left_L/units::cm+1e-9)/50e3) << std::endl;
	return false;
      }
    }else{
      bool flag_fix_end = false;
      if (exit_L < 35*units::cm || (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0&& left_L > 2*units::cm) flag_fix_end = true;
      
      if (left_L < 8*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3)< 1.5 ||
	  left_L < 6*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.7 ||
	  left_L < 5*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.8 ||
	  left_L < 3*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.9){
	left_L = 0;
	kink_num = dQ.size();
	exit_L = 40*units::cm;
	flag_fix_end = false;
      }
    
     
      bool flag_pass = false;

      // std::cout << flag_other_clusters << std::endl;
      
      if (!flag_other_clusters){
	if (left_L < 40*units::cm) {
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
	
	if (left_L < 20*units::cm){
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
	    }
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm , 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm , 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
      }else{
	if (flag_fix_end)
	  flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm, true);
	else
	  flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm, 0., 35*units::cm, true);
      }
      
      if (flag_pass) {
	main_cluster->clear_fit_tracks();
	main_cluster->search_other_tracks(ct_point_cloud, global_wc_map, flash_time*units::microsecond);

	if (check_other_tracks(main_cluster, offset_x)){
	  std::cout << "Mid Point Tracks" << std::endl;
	  return false;
	}
	
	//	std::cout << main_cluster->get_fit_tracks().size() << std::endl;
	
	if (!detect_proton(main_cluster, kink_num)) return true;
      }
    }
  }
  
  // backward check ...
  if (flag_double_end){
    std::cout << "Backward check! " << std::endl;
    // regular crawling ...
    main_cluster->do_rough_path( last_wcp, first_wcp);
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond, false);
    Point mid_p = main_cluster->adjust_rough_path(); 
    // fitting trajectory and dQ/dx...
    main_cluster->collect_charge_trajectory(ct_point_cloud); 
    main_cluster->do_tracking(ct_point_cloud, global_wc_map, flash_time*units::microsecond);

    // check both end points for TGM ...
    WCP::PointVector& pts = main_cluster->get_fine_tracking_path();
    std::vector<double>& dQ = main_cluster->get_dQ();
    std::vector<double>& dx = main_cluster->get_dx();
    
    int kink_num = find_first_kink(main_cluster);
    
    double left_L = 0;
    double left_Q = 0;
    double exit_L = 0;
    double exit_Q = 0;
    for (size_t i=0;i!=kink_num;i++){
      exit_L += dx.at(i);
      exit_Q += dQ.at(i);
    }
    for (size_t i = kink_num; i!=dx.size(); i++){
      left_L += dx.at(i);
      left_Q += dQ.at(i);
    }
    
    std::cout << "Left: " << exit_L/units::cm << " " << left_L/units::cm << " " << (left_Q/(left_L/units::cm+1e-9))/50e3 << " " << (exit_Q/(exit_L/units::cm+1e-9)/50e3) << std::endl;
    
    if ( (!inside_fiducial_volume(pts.front(),offset_x)) && (!inside_fiducial_volume(pts.back(),offset_x))){
      bool flag_TGM_anode = false;
      
      if ((pts.back().x < 2*units::cm || pts.front().x < 2*units::cm) && kink_num >=0 && kink_num < pts.size()){ // at Anode ...
	if (pts.at(kink_num).x < 6*units::cm){
	  TVector3 v10(pts.back().x-pts.at(kink_num).x,pts.back().y-pts.at(kink_num).y,pts.back().z-pts.at(kink_num).z);
	  TVector3 v20(pts.front().x-pts.at(kink_num).x,pts.front().y-pts.at(kink_num).y,pts.front().z-pts.at(kink_num).z);
	  if (fabs(v10.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v10.Mag()>15*units::cm || fabs(v20.Angle(drift_dir)/3.1415926*180.-90)<12.5 && v20.Mag()>15*units::cm)
	    flag_TGM_anode = true;
	  //std::cout << v10.Angle(drift_dir)/3.1415926*180. << " " << v20.Angle(drift_dir)/3.1415926*180. << std::endl;
	}
      }
      if ((exit_L < 3*units::cm || left_L < 3*units::cm) || flag_TGM_anode){
	std::cout << "TGM: " << pts.front() << " " << pts.back() << std::endl;
	//event_type |= 1UL << 3;
	return false;
      }
    }
  
    if (left_L > 40*units::cm || left_L > 7.5*units::cm && (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0){
      std::cout << "Mid Point A " << inside_dead_region(mid_p) << " " << mid_p << " " << left_L  << " " << (left_Q/(left_L/units::cm+1e-9)/50e3) << std::endl;
      return false;
    }else{
      bool flag_fix_end = false;
      if (exit_L < 35*units::cm || (left_Q/(left_L/units::cm+1e-9))/50e3 > 2.0 && left_L > 2*units::cm) flag_fix_end = true;
      
      if (left_L < 8*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3)< 1.5 ||
	  left_L < 6*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.7 ||
	  left_L < 3*units::cm && (left_Q/(left_L/units::cm+1e-9)/50e3) < 1.9){
	left_L = 0;
	kink_num = dQ.size();
	exit_L = 40*units::cm;
	flag_fix_end = false;
      }

    
    
      bool flag_pass = false;
      if (!flag_other_clusters){
	if (left_L < 40*units::cm) {
	  if (flag_fix_end){
	    flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	  }else{
	    flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 35*units::cm) ||
	      eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 35*units::cm);
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 40*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
	
	if (left_L < 20*units::cm){
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm, 3.*units::cm, 35*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 35*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 35*units::cm);
	    }
	  }
	  
	  if (!flag_pass){
	    if (flag_fix_end){
	      flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm , 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 5*units::cm , 3.*units::cm, 15*units::cm);
	    }else{
	      flag_pass = eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 0., 15*units::cm) ||
		eval_stm(main_cluster, kink_num, 20*units::cm - left_L, 3.*units::cm, 15*units::cm);
	    }
	  }
	}
      }else{
	if (flag_fix_end)
	  flag_pass = eval_stm(main_cluster, kink_num, 5*units::cm, 0., 35*units::cm, true);
	else
	  flag_pass = eval_stm(main_cluster, kink_num, 40*units::cm, 0., 35*units::cm, true);
      }

      if (flag_pass) {
	main_cluster->clear_fit_tracks();
	main_cluster->search_other_tracks(ct_point_cloud, global_wc_map, flash_time*units::microsecond);
	if (check_other_tracks(main_cluster, offset_x)){
	  std::cout << "Mid Point Tracks" << std::endl;
	  return false;
	}
	
	if (!detect_proton(main_cluster, kink_num)) return true;
      }
    }
  }
  
  // check 5512-209-10491
  // if (inside_dead_region(mid_p)) return true;
  //  

  // end check ...
  std::cout << "Mid Point " << std::endl;  
  return false;
}

