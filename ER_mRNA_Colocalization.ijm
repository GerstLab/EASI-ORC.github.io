cell_dir = getDirectory("Cell mask images directory:");
cell_list = getFileList(cell_dir);
DAPI_dir = getDirectory("DAPI mask images directory:");
DAPI_list = getFileList(DAPI_dir);
smFISH_dir = getDirectory("smFISH images directory (probe channel):");
smFISH_list = getFileList(smFISH_dir);
nER_dir = getDirectory("Perinuclear ER mask images directory:");
nER_list = getFileList(nER_dir);
cER_dir = getDirectory("Cortical ER mask images directory:");
cER_list = getFileList(cER_dir);
results_path = getDirectory("Results Tables Location:");
prominence_val = parseInt(getString("Prominance value for the smFISH signal:", "100"));
smFISH_signal_cutoff = parseInt(getString("Approximate radius of an smFISH signal (in pixels):", "2"));
ER_signal_min = parseInt(getString("Minimum amount of ER signal within a cell:", "8"));
ER_signal_max = parseInt(getString("Maximum amount of ER signal within a cell:", "45"));


//setBatchMode(true);
for (img = 0; img < cell_list.length; img++) {
	//Make sure file is not a folder (won't go into subfolders)
	cell_img_path = cell_dir + cell_list[img];
	dapi_img_path = DAPI_dir + DAPI_list[img];
	smFISH_img_path = smFISH_dir + smFISH_list[img];
	nER_img_path = nER_dir + nER_list[img];
	cER_img_path = cER_dir + cER_list[img];
	if (File.isFile(cell_img_path)) {
		open(cell_img_path);
		//Create array to mark each pixel in image for mRNA maxima point
		img_width = getWidth();
		img_height = getHeight();
		maxima_matrix = newArray(img_width * img_height);
		//Dilate to improve cER cover and use Watershed to seperate close cells
		setOption("BlackBackground", false);
		run("Dilate", "stack");
		run("Watershed", "stack");
		//Create folders for seperated slices and maps
		cell_seperated_folder = cell_dir + "seperated_planes";
		File.makeDirectory(cell_seperated_folder);
		//Go over cells image and create a file from each slice in a new directory
	 	run("Image Sequence... ", "dir=[" + cell_seperated_folder + "] format=TIFF digits=2");
		run("Close All");
		// Go over dapi image and create a file from each slice in a new directory
		if (File.isFile(dapi_img_path)) {
			//Create folders for seperated slices
			dapi_seperated_folder = DAPI_dir + "seperated_planes";
			File.makeDirectory(dapi_seperated_folder);
	 		open(dapi_img_path);
	 		run("Image Sequence... ", "dir=[" + dapi_seperated_folder + "] format=TIFF digits=2");
			run("Close All");
		}
		// Go over smFISH image and create a file from each slice in a new directory
		if (File.isFile(smFISH_img_path)) {
			//Create folders for seperated slices
			smFISH_seperated_folder = smFISH_dir + "seperated_planes";
			File.makeDirectory(smFISH_seperated_folder);
	 		open(smFISH_img_path);
	 		run("Image Sequence... ", "dir=[" + smFISH_seperated_folder + "] format=TIFF digits=2");
			run("Close All");
		}
		// Go over nER image and create a file from each slice in a new directory
		if (File.isFile(nER_img_path)) {
			//Create folders for seperated slices
			nER_seperated_folder = nER_dir + "seperated_planes";
			File.makeDirectory(nER_seperated_folder);
	 		open(nER_img_path);
	 		run("Image Sequence... ", "dir=[" + nER_seperated_folder + "] format=TIFF digits=2");
			run("Close All");
		}
		// Go over cER image and create a file from each slice in a new directory
		if (File.isFile(cER_img_path)) {
			//Create folders for seperated slices
			cER_seperated_folder = cER_dir + "seperated_planes";
			File.makeDirectory(cER_seperated_folder);
	 		open(cER_img_path);
	 		run("Image Sequence... ", "dir=[" + cER_seperated_folder + "] format=TIFF digits=2");
			run("Close All");
		}
		//Go over seperate slices and create ROIs from the cells image
		cell_plane_list = getFileList(cell_seperated_folder);
		dapi_plane_list = getFileList(dapi_seperated_folder);
		smFISH_plane_list = getFileList(smFISH_seperated_folder);
		nER_plane_list = getFileList(nER_seperated_folder);
		cER_plane_list = getFileList(cER_seperated_folder);
		cell_slice = cell_seperated_folder + File.separator + cell_plane_list[0];
		open(cell_slice);
		run("Analyze Particles...", "exclude add");
		//Get total number of cells for image and create arrays for data points
		run("Close All");
		nCells = roiManager("count");
	 	tot_mRNA_per_cell = newArray(nCells);
	 	tot_colo_nER = newArray(nCells);
	 	tot_colo_cER = newArray(nCells);
	 	tot_no_colo = newArray(nCells);
	 	ER_signal_size = newArray(nCells);
		for (plane = 0; plane < cell_plane_list.length; plane++) {
			dapi_signal_list = newArray(nCells);
			edge_cells_list = newArray(nCells);
			dapi_slice = dapi_seperated_folder + File.separator + dapi_plane_list[plane];
			smFISH_slice = smFISH_seperated_folder + File.separator + smFISH_plane_list[plane];
			nER_slice = nER_seperated_folder + File.separator + nER_plane_list[plane];
			cER_slice = cER_seperated_folder + File.separator + cER_plane_list[plane];
			open(cER_slice);
			current_cER = getTitle();
			open(nER_slice);
			current_nER = getTitle();
			if (isOpen("ROI Manager")) {
				if (nCells == 0) {
					continue;
				}
				roiManager("select", Array.getSequence(roiManager("count")));
				roiManager("delete");
			}
			open(cell_slice);
			//Exclude any ROI whithin edges (2% of pixels in img)
			max_x = img_width * 0.98;
			min_x = img_width * 0.02;
			max_y = img_height * 0.98;
			min_y = img_height * 0.02;
			run("Analyze Particles...", "exclude add");
			filter_list = newArray();
			for (roi = 0; roi < nCells; roi++) {
				roi_xpoints = newArray();
				roi_ypoints = newArray();
				roiManager("select", roi);
				Roi.getContainedPoints(roi_xpoints, roi_ypoints);
				for (pix = 0; pix < roi_xpoints.length; pix++) {
					if (roi_xpoints[pix] > max_x || roi_xpoints[pix] < min_x || roi_ypoints[pix] > max_y || roi_ypoints[pix] < min_y) {
						edge_cells_list[roi] = 1;
						break;
					}
				}
			}
			//Register DAPI signal in each cell (disregard planes where there's no DAPI signal in cell)
			open(dapi_slice);
			run("Set Measurements...", "area_fraction redirect=None decimal=5");
			setThreshold(255, 255, "raw");
			roiManager("select", Array.getSequence(nCells));
			roiManager("measure");
			cells_table = "Cells Table "+ plane;
			IJ.renameResults(cells_table);
			//Go over each result row and and mark any cells with dapi signal
			dapi_cutoff = 0.15;
			for (roi = 0; roi < nCells; roi++) {
				if (Table.get("%Area", roi, cells_table) < dapi_cutoff) {
					dapi_signal_list[roi] = 1;
				}
			}
			//Register % of ER signal within cell (disregard planes where signal is too great/small)
			selectWindow(current_cER);
			roiManager("select", Array.getSequence(nCells));
			roiManager("measure");
			cER_signal_table = "cER Signal " + plane;
			IJ.renameResults(cER_signal_table);
			selectWindow(current_nER);
			roiManager("measure");			
			nER_signal_table = "nER Signal " + plane;
			IJ.renameResults(nER_signal_table);
			for (roi = 0; roi < nCells; roi++) {
				if (Table.get("%Area", roi, cER_signal_table) +  Table.get("%Area", roi, nER_signal_table) < ER_signal_min ||Table.get("%Area", roi, cER_signal_table) +  Table.get("%Area", roi, nER_signal_table) > ER_signal_max) {
					ER_signal_size[roi] = 1;
				}
			}
			//At this point, we know which of the cells are at the dges, have a dapi signal, and the % of ER signals in each cell
			//Find Maxima for smFISH signals and check in which cell they are located
			open(smFISH_slice);
			roiManager("select", Array.getSequence(roiManager("count")));
			roiManager("Combine");
			wait(500);
			run("Clear Results");
			setBackgroundColor(0, 0, 0);
			run("Clear Outside");
			run("Find Maxima...", "prominence=" + prominence_val + " strict output=List");
			mRNA_table = "mRNA Table " + plane;
			setResult("Duplicate", 0, 0);
			IJ.renameResults(mRNA_table);
			num_of_mRNA = Table.size(mRNA_table);
			for (mRNA = 0; mRNA < num_of_mRNA; mRNA++) {
				for (cell = 0; cell < nCells; cell++) {
					if (dapi_signal_list[cell] == 1 || edge_cells_list[cell] == 1 || ER_signal_size[cell] == 1) {
						continue;
					}
					roiManager("select", cell);
					current_X = Table.get("X", mRNA, mRNA_table);
					current_Y = Table.get("Y", mRNA, mRNA_table);
					if (Roi.contains(current_X, current_Y)) {
						//Create lists of pixels in bounding rectangle of circle around point (to check duplicates)
						rect_x = current_X - smFISH_signal_cutoff;
						rect_y = current_Y - smFISH_signal_cutoff;
						rect_size = (smFISH_signal_cutoff * 2);
						radius_sqr = smFISH_signal_cutoff * smFISH_signal_cutoff;
						rect_xs = newArray(rect_size * rect_size);
						rect_ys = newArray(rect_size * rect_size);
						for (x = 0; x < rect_size; x++) {
							for (y = 0; y < rect_size; y++) {
								rect_xs[y * rect_size + x] = rect_x + x;
								rect_ys[y * rect_size + x] = rect_y + y;
							}
						}
						//Calc distance of each point in bounding sqaure from maxima and add points within circle with radius of cutoff to array
						near_points_list = newArray();
						for (point = 0; point < rect_xs.length; point++) {
							distance_x = rect_xs[point] - current_X;
							distance_y = rect_ys[point] - current_Y;
							distance_sqr = (distance_x * distance_x) + (distance_y * distance_y);
							if (distance_sqr < radius_sqr) {
								near_point = newArray();
								pixel_location = (rect_ys[point] * img_width) + rect_xs[point];
								near_point[0] = pixel_location;
								near_points_list = Array.concat(near_points_list, near_point);
							}
						}
						//Check each near point, to make sure it wasn't identified before
						for (point = 0; point < near_points_list.length; point++) {
							if (maxima_matrix[near_points_list[point]] == 0 && point == near_points_list.length - 1) {
								maxima_matrix[(current_Y * img_width) + current_X] += 1;
								Table.set("Cell In Slice " + plane, mRNA, cell, mRNA_table);
								tot_mRNA_per_cell[cell]++;
								//After finding and saving the cell in which the mRNA is, we look at colocalization with c/nER
								//Use ability to draw selection circles on ER mask images
								selectWindow(current_cER);
								//Create circle selection with radius of cutoff given by user, around maxima pixel
								makeOval(current_X - smFISH_signal_cutoff, current_Y - smFISH_signal_cutoff, smFISH_signal_cutoff * 2, smFISH_signal_cutoff * 2);
								run("Clear Results");
								run("Set Measurements...", "mean redirect=None decimal=2");
								run("Measure");
								selectWindow(current_nER);
								//Create circle selection with radius of cutoff given by user, around maxima pixel
								makeOval(current_X - smFISH_signal_cutoff, current_Y - smFISH_signal_cutoff, smFISH_signal_cutoff * 2, smFISH_signal_cutoff * 2);
								run("Measure");
									if (getResult("Mean", 0) != 0 && getResult("Mean", 1) != 0) {
										if (getResult("Mean", 0) > getResult("Mean", 1)) {
											Table.set("cER_Colo_Slice" + plane, mRNA, cell, mRNA_table);
											tot_colo_cER[cell]++;
											continue;											
										}
										else if (getResult("Mean", 0) <= getResult("Mean", 1)) {
											Table.set("nER_Colo_Slice" + plane, mRNA, cell, mRNA_table);
											tot_colo_nER[cell]++;
											continue;
										}
									}
									else if (getResult("Mean", 0) != 0) {
										Table.set("cER_Colo_Slice" + plane, mRNA, cell, mRNA_table);
										tot_colo_cER[cell]++;
										continue;										
									}
									else if (getResult("Mean", 1) != 0) {
										Table.set("nER_Colo_Slice" + plane, mRNA, cell, mRNA_table);
										tot_colo_nER[cell]++;
										continue;
									}
									else {
										Table.set("No_Colo_Slice" + plane, mRNA, cell, mRNA_table);
										tot_no_colo[cell]++;
										continue;
									}
							}
							else if (maxima_matrix[near_points_list[point]] > 0) {
							Table.set("Duplicate", mRNA, 1, mRNA_table);
							break;
							}
						}
						//Break Roi loop (mRNA already catagorized) and go to next mRNA
						break;
					}
					//Break cell loop (mRNA already catagorized) and go to next mRNA
					//End of Roi loop
				}
				//End of cell loop
			}
			//End of mRNA loop
			close("*");
			close("Roi Manager");
			close(mRNA_table);
			close(cells_table);
			close(cER_signal_table);
			close(nER_signal_table);
		}
		//End of Plane loop
		close("*");
		
		
		//Remove single plane files
		//Cells
		for (sep = 0; sep < cell_plane_list.length; sep++) {
			img_to_del = cell_seperated_folder + File.separator + cell_plane_list[sep];
			File.delete(img_to_del);
		}
		File.delete(cell_seperated_folder);
		//DAPI
		for (sep = 0; sep < dapi_plane_list.length; sep++) {
			img_to_del = dapi_seperated_folder + File.separator + dapi_plane_list[sep];
			File.delete(img_to_del);
		}
		File.delete(dapi_seperated_folder);
		//smFISH
		for (sep = 0; sep < smFISH_plane_list.length; sep++) {
			img_to_del = smFISH_seperated_folder + File.separator + smFISH_plane_list[sep];
			File.delete(img_to_del);
		}
		File.delete(smFISH_seperated_folder);
		//nER
		for (sep = 0; sep < nER_plane_list.length; sep++) {
			img_to_del = nER_seperated_folder + File.separator + nER_plane_list[sep];
			File.delete(img_to_del);
		}
		File.delete(nER_seperated_folder);
		//cER
		for (sep = 0; sep < cER_plane_list.length; sep++) {
			img_to_del = cER_seperated_folder + File.separator + cER_plane_list[sep];
			File.delete(img_to_del);
		}
		File.delete(cER_seperated_folder);
	}
	//Sum and save data for img
	results_table_name = cell_list[img] + " ER Colocalization.csv";
	Table.create(results_table_name);
	Table.setColumn("Cell #", Array.getSequence(tot_mRNA_per_cell.length), results_table_name);
	Table.setColumn("Total mRNA per Cell", tot_mRNA_per_cell, results_table_name);
	Table.setColumn("Total Colocolized With nER", tot_colo_nER, results_table_name);
	Table.setColumn("Total Colocolized With cER", tot_colo_cER, results_table_name);
	Table.setColumn("Total Not Colocolized with ER", tot_no_colo, results_table_name);
	//Add edge filtering column and remove appropriate lines from results table
	Table.setColumn("Edge Cells", edge_cells_list, results_table_name);
	for (row = nCells - 1; row >= 0; row--) {
		if (Table.get("Edge Cells", row, results_table_name) == 1 ) {
			Table.deleteRows(row, row, results_table_name);
		}
	}
	Table.deleteColumn("Edge Cells", results_table_name);
	selectWindow(results_table_name);
	//Save results as csv
	saveAs("Results", results_path + results_table_name);
	run("Close All");
	close(results_table_name);
}
//End of img loop
print ("All images done. Results saved under " + results_path);
