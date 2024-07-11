//get the location of cells masks images
Dialog.create("Cells masks location");
Dialog.addMessage("Please input the cells mask images directory:");
Dialog.show();
cell_dir = getDirectory("Cells mask images directory:");
cell_list = getFileList(cell_dir);
//get the location of smFISH csv files
Dialog.create("smFISH csv files location");
Dialog.addMessage("Please input the smFISH spots csv files directory:");
Dialog.show();
smFISH_dir = getDirectory("smFISH csv files directory (created by RSFISH):");
smFISH_list = getFileList(smFISH_dir);
//ask user if ER is segmented and get proper mask images
Dialog.create("ER Sub-segmentation");
Dialog.addChoice("Was the organelle sub-segmented (Select 'No' if you have one organelle mask)?", newArray("Yes", "No"));
Dialog.show();
ER_seg_check = Dialog.getChoice();
if (ER_seg_check == "Yes") {
	Dialog.create("Near masks location");
	Dialog.addMessage("Please input the near organelle mask images directory:");
	Dialog.show();
	nER_dir = getDirectory("Near organelle mask images directory:");
	nER_list = getFileList(nER_dir);
	Dialog.create("Far masks location");
	Dialog.addMessage("Please input the far organelle mask images directory:");
	Dialog.show();
	cER_dir = getDirectory("Far mask images directory:");
	cER_list = getFileList(cER_dir);
	//get the location of proximity marker mask images 
	Dialog.create("Proximity marker masks location");
	Dialog.addMessage("Please input the proximity marker mask images directory:");
	Dialog.show();
	DAPI_dir = getDirectory("Proximity marker mask images directory:");
	DAPI_list = getFileList(DAPI_dir);
	prox_mark_filter = parseFloat(getString("Choose a minimum threshold for proximity marker cell coverage per plane", "0.05"));
}
if (ER_seg_check == "No") {
	Dialog.create("Organelle masks location");
	Dialog.addMessage("Please input the Organelle mask images directory:");
	Dialog.show();
	whole_ER_dir = getDirectory("Organelle mask images directory:");
	whole_ER_list = getFileList(whole_ER_dir);
}
smFISH_signal_cutoff = parseFloat(getString("Diameter (in pixels) for spot colocaliztion calculation", "1"));
if (smFISH_signal_cutoff <= 0) {
	smFISH_signal_cutoff = parseFloat(getString("Only values greater than 0 are possible ", "1"));
}

Dialog.create("Results tables location");
Dialog.addMessage("Please input the directory you wish results tables will be saved:");
Dialog.show();
results_path = getDirectory("Results Tables Location:");


for (img = 0; img < cell_list.length; img++) {
	//Make sure file is not a folder (won't go into subfolders)
	cell_img_path = cell_dir + cell_list[img];
	if (ER_seg_check == "Yes") {
		dapi_img_path = DAPI_dir + DAPI_list[img];
		nER_img_path = nER_dir + nER_list[img];
		cER_img_path = cER_dir + cER_list[img];
	}
	if (ER_seg_check == "No") {
		ER_img_path = whole_ER_dir + whole_ER_list[img];
	}
	smFISH_table = smFISH_dir + smFISH_list[img];
	smFISH_counter = 0;
	if (File.isFile(cell_img_path)) {
		open(smFISH_table);
		Table.sort("z");
		open(cell_img_path);
		ncross_section = nSlices;
		//Dilate to improve cellular periphery cover and use Watershed to seperate close cells
		setOption("BlackBackground", false);
		run("Dilate", "stack");
		run("Watershed", "stack");
		//Create folders for seperated slices and maps
		cell_seperated_folder = cell_dir + "seperated_planes";
		File.makeDirectory(cell_seperated_folder);
		//Go over cells image and create a file from each slice in a new directory
	 	run("Image Sequence... ", "dir=[" + cell_seperated_folder + "] format=TIFF digits=2");
		run("Close All");
		if (ER_seg_check == "Yes") {
			// Go over dapi image and create a file from each slice in a new directory
			if (File.isFile(dapi_img_path)) {
				//Create folders for seperated slices
				dapi_seperated_folder = DAPI_dir + "seperated_planes";
				File.makeDirectory(dapi_seperated_folder);
	 			open(dapi_img_path);
	 			run("Image Sequence... ", "dir=[" + dapi_seperated_folder + "] format=TIFF digits=2");
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
		}
		// Go over non-segmented ER/other organelle images and create a file from each slice in a new directory
		if (ER_seg_check == "No") {
			ER_seperated_folder = whole_ER_dir + "seperated_planes";
			File.makeDirectory(ER_seperated_folder);
	 		open(ER_img_path);
	 		run("Image Sequence... ", "dir=[" + ER_seperated_folder + "] format=TIFF digits=2");
			run("Close All");
		}

		//Go over seperate slices and create ROIs from the cells image
		cell_plane_list = getFileList(cell_seperated_folder);
		if (ER_seg_check == "Yes") {
			dapi_plane_list = getFileList(dapi_seperated_folder);
			nER_plane_list = getFileList(nER_seperated_folder);
			cER_plane_list = getFileList(cER_seperated_folder);
		}
		if (ER_seg_check == "No") {
			ER_plane_list = getFileList(ER_seperated_folder);
		}
		cell_slice = cell_seperated_folder + File.separator + cell_plane_list[0];
		open(cell_slice);
		img_width = getWidth();
		img_height = getHeight();
		run("Analyze Particles...", "exclude add");
		//Get total number of cells for image and create arrays for mRNA spot counting of each cell
		run("Close All");
		nCells = roiManager("count");
	 	tot_mRNA_per_cell = newArray(nCells);
	 	spot_data_per_cell = newArray(nCells);
	 	if (ER_seg_check == "Yes") {
	 		tot_colo_nER = newArray(nCells);
	 		tot_colo_cER = newArray(nCells);
	 	}
	 	if (ER_seg_check == "No") {
	 		tot_colo_ER = newArray(nCells);
	 		zero_column = newArray(nCells);
	 	}
	 	tot_no_colo = newArray(nCells);
	 	//Create table to save  ER signal values for each plane
	 	ER_signal_size_name = "ER Signal Size Per Plane";
		Table.create(ER_signal_size_name);
	 	}
		for (plane = 0; plane < cell_plane_list.length; plane++) {
			edge_cells_list = newArray(nCells);
			if (ER_seg_check == "Yes") {
				dapi_signal_list = newArray(nCells);
				dapi_slice = dapi_seperated_folder + File.separator + dapi_plane_list[plane];
				nER_slice = nER_seperated_folder + File.separator + nER_plane_list[plane];
				cER_slice = cER_seperated_folder + File.separator + cER_plane_list[plane];
				open(cER_slice);
				current_cER = getTitle();
				open(nER_slice);
				current_nER = getTitle();
			}
			if (ER_seg_check == "No") {
				ER_slice = ER_seperated_folder + File.separator + ER_plane_list[plane];
				open(ER_slice);
				current_ER = getTitle();
			}
			if (isOpen("ROI Manager")) {
				if (nCells == 0) {
					continue;
				}
				roiManager("select", Array.getSequence(roiManager("count")));
				roiManager("delete");
			}
			open(cell_slice);
			//Exclude any cell ROI whithin edges (2% of pixels in img)
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
			run("Set Measurements...", "area_fraction redirect=None decimal=5");
			//Array for saving current plane's ER signal %Area in each cell
			current_ER_signal_size = newArray(nCells);
			if (ER_seg_check == "Yes") {
				//Register DAPI signal in each cell (disregard planes where there's no DAPI signal in cell - only for sub-segmented ER)
				open(dapi_slice);
				setThreshold(255, 255, "raw");
				roiManager("select", Array.getSequence(nCells));
				roiManager("measure");
				cells_table = "Cells Table "+ plane;
				IJ.renameResults(cells_table);
				//Go over each result row and and mark any cells with dapi signal
				for (roi = 0; roi < nCells; roi++) {
					if (Table.get("%Area", roi, cells_table) < prox_mark_filter) {
						dapi_signal_list[roi] = 1;
					}
				}
				//Register % of ER signal within cell (will allow to filter planes with too little or to mich ER)
				selectWindow(current_cER);
				roiManager("select", Array.getSequence(nCells));
				roiManager("measure");
				cER_signal_table = "cER Signal " + plane;
				IJ.renameResults(cER_signal_table);
				selectWindow(current_nER);
				roiManager("measure");
				nER_signal_table = "nER Signal " + plane;
				IJ.renameResults(nER_signal_table);
				for (result = 0; result < nCells; result++) {
					current_ER_signal_size[result] = getResult("%Area", result, nER_signal_table) + getResult("%Area", result, cER_signal_table);
				}
			}
			if (ER_seg_check == "No") {
				//Register % of organelle signal within cell (disregard planes where signal is too great/small)
				selectWindow(current_ER);
				roiManager("select", Array.getSequence(nCells));
				roiManager("measure");
				ER_signal_table = "ER Signal" + plane;
				IJ.renameResults(ER_signal_table);
				for (result = 0; result < nCells; result++) {
					current_ER_signal_size[result] = getResult("%Area", result, ER_signal_table);
				}
			}
			//Save ER size values for this plane
			Table.setColumn("Organelle Signal Size Cross Section " + (plane + 1), current_ER_signal_size, ER_signal_size_name);
			//Go over mRNAs in this slice, find each's cell and measure colocalization
			num_of_mRNA = Table.size(smFISH_list[img]);
			for (mRNA = smFISH_counter; mRNA < num_of_mRNA; mRNA++) {
				//Check "z" column for spot's slice. If z is smaller than current slice + 1, the spot is in current slice.
				if (round(Table.get("z", mRNA, smFISH_list[img])) > (plane + 1)) {
					break;
				}
				smFISH_counter += 1;
				//Skip smFISH spots with 'z' coordiantes smaller than 0 (no organelle data exsists for them)
				if (Table.get("z", mRNA, smFISH_list[img]) < 0 || Table.get("z", mRNA, smFISH_list[img]) > (ncross_section - 1)) {
					continue;
				}
				//Check for ER, DAPI and edge limits for each cell
				for (cell = 0; cell < nCells; cell++) {
					if (edge_cells_list[cell] == 1) {
						continue;
					}
					if (ER_seg_check == "Yes") {
						if (dapi_signal_list[cell] == 1) {
							continue;
						}
					}
					roiManager("select", cell);
					current_X = Table.get("x", mRNA, smFISH_list[img]);
					current_Y = Table.get("y", mRNA, smFISH_list[img]);
					current_Z = round(Table.get("z", mRNA, smFISH_list[img]));
					current_intesity = Table.get("intensity", mRNA, smFISH_list[img]);
					if (Roi.contains(current_X, current_Y)) {
						tot_mRNA_per_cell[cell]++;
						//Use selection circles on organelle mask image - sub-segmented ER
						if (ER_seg_check == "Yes") {
							selectWindow(current_cER);
							//Create circle selection with diameter of cutoff given by user, around smFISH spot pixel
							makeOval(current_X - smFISH_signal_cutoff/2 + 0.5, current_Y - smFISH_signal_cutoff/2 + 0.5, smFISH_signal_cutoff, smFISH_signal_cutoff);
							run("Clear Results");
							run("Set Measurements...", "mean redirect=None decimal=2");
							run("Measure");
							selectWindow(current_nER);
							//Create circle selection with diameter of cutoff given by user, around smFISH spot pixel
							makeOval(current_X - smFISH_signal_cutoff/2 + 0.5, current_Y - smFISH_signal_cutoff/2 + 0.5, smFISH_signal_cutoff, smFISH_signal_cutoff);
							run("Measure");
							//Check if there's ER signal near the smFISH spot and save spot data
							if (getResult("Mean", 0) != 0 && getResult("Mean", 1) != 0) {
								//Case where signals of two organelles/sub-organelles are equal
								equal_radius = 0.5;
								while (getResult("Mean", 0) == getResult("Mean", 1)) {
									equal_radius += 0.5;
									selectWindow(current_cER);
									//Create circle selection with diameter of 1 and increasing, around smFISH spot pixel
									makeOval(current_X - equal_radius/2 + 0.5, current_Y - equal_radius/2 + 0.5, equal_radius, equal_radius);
									run("Clear Results");
									run("Set Measurements...", "mean redirect=None decimal=2");
									run("Measure");
									selectWindow(current_nER);
									//Create circle selection with diameter of 1 and increasing, around smFISH spot pixel
									makeOval(current_X - equal_radius/2 + 0.5, current_Y - equal_radius/2 + 0.5, equal_radius, equal_radius);
									run("Measure");
									continue;
								}
								if (getResult("Mean", 0) > getResult("Mean", 1)) {
									tot_colo_cER[cell]++;
									if (spot_data_per_cell[cell] == 0) {
										spot_data_per_cell[cell] = "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_far],";
									}
									else {
										spot_data_per_cell[cell] = spot_data_per_cell[cell] + "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_far],";
									}
									continue;
								}
								else if (getResult("Mean", 0) < getResult("Mean", 1)) {
									tot_colo_nER[cell]++;
									if (spot_data_per_cell[cell] == 0) {
										spot_data_per_cell[cell] = "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_near],";
									}
									else {
										spot_data_per_cell[cell] = spot_data_per_cell[cell] + "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_near],";
									}
									continue;
								}
							}
							else if (getResult("Mean", 0) != 0) {
							tot_colo_cER[cell]++;
								if (spot_data_per_cell[cell] == 0) {
									spot_data_per_cell[cell] = "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_far],";
								}
								else {
									spot_data_per_cell[cell] = spot_data_per_cell[cell] + "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_far],";
								}
							continue;
							}
							else if (getResult("Mean", 1) != 0) {
								tot_colo_nER[cell]++;
								if (spot_data_per_cell[cell] == 0) {
									spot_data_per_cell[cell] = "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_near],";
								}
								else {
									spot_data_per_cell[cell] = spot_data_per_cell[cell] + "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_near],";
								}
								continue;
							}
							else {
								tot_no_colo[cell]++;
								if (spot_data_per_cell[cell] == 0) {
									spot_data_per_cell[cell] = "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",nc],";
								}
								else {
									spot_data_per_cell[cell] = spot_data_per_cell[cell] + "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",nc],";
								}
								continue;
							}
						}
						//Use ability to draw selection circles on organelle mask image - single mask organelle
						if (ER_seg_check == "No") {
							zero_column[cell] = 0;
							selectWindow(current_ER);
							//Create circle selection with diameter of cutoff given by user, around smFISH spot pixel
							makeOval(current_X - smFISH_signal_cutoff/2 + 0.5, current_Y - smFISH_signal_cutoff/2 + 0.5, smFISH_signal_cutoff, smFISH_signal_cutoff);
							run("Clear Results");
							run("Set Measurements...", "mean redirect=None decimal=2");
							run("Measure");
							//Check if there's ER signal near the smFISH spot
							if (getResult("Mean", 0) != 0) {
								tot_colo_ER[cell]++;
								if (spot_data_per_cell[cell] == 0) {
								spot_data_per_cell[cell] = "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_near],";
								}
								else {
									spot_data_per_cell[cell] = spot_data_per_cell[cell] + "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",organelle_near],";
								}
								continue;
							}
							else {
								tot_no_colo[cell]++;
								if (spot_data_per_cell[cell] == 0) {
								spot_data_per_cell[cell] = "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",nc],";
								}
								else {
									spot_data_per_cell[cell] = spot_data_per_cell[cell] + "[" + current_X + "," + current_Y + "," + current_Z + "," + current_intesity + ",nc],";
								}
								continue;
							}
						}

					//End of Roi loop
					}
				//End of cell loop
				}
			//End of mRNA loop
			}
			close("*");
			close("Roi Manager");
			if (ER_seg_check == "Yes") {
				close(cells_table);
				close(cER_signal_table);
				close(nER_signal_table);
			}
			if (ER_seg_check == "No") {
				close(ER_signal_table);
			}
		//End of Plane loop
		}
		close("*");
		//Remove single plane files
		//Cells
		for (sep = 0; sep < cell_plane_list.length; sep++) {
			img_to_del = cell_seperated_folder + File.separator + cell_plane_list[sep];
			File.delete(img_to_del);
		}
		File.delete(cell_seperated_folder);
		if (ER_seg_check == "Yes") {
			//DAPI
			for (sep = 0; sep < dapi_plane_list.length; sep++) {
				img_to_del = dapi_seperated_folder + File.separator + dapi_plane_list[sep];
				File.delete(img_to_del);
			}
			File.delete(dapi_seperated_folder);
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
		if (ER_seg_check == "No") {
			for (sep = 0; sep < ER_plane_list.length; sep++) {
				img_to_del = ER_seperated_folder + File.separator + ER_plane_list[sep];
				File.delete(img_to_del);
			}
			File.delete(ER_seperated_folder);
		}
	//Sum and save data for img
	//create the final results table
	results_table_name = cell_list[img] + "Colocalization.csv";
	Table.create(results_table_name);
	//Add cells index column
	Table.setColumn("Cell #", Array.getSequence(tot_mRNA_per_cell.length), results_table_name);
	//Add number of mRNA and colocalization columns
	Table.setColumn("Total mRNA per Cell", tot_mRNA_per_cell, results_table_name);
	if (ER_seg_check == "Yes") {
		Table.setColumn("Total Colocolized With Organelle Near", tot_colo_nER, results_table_name);
		Table.setColumn("Total Colocolized With Organelle Far", tot_colo_cER, results_table_name);
	}
	if (ER_seg_check == "No") {
		Table.setColumn("Total Colocolized With Organelle Near", tot_colo_ER, results_table_name);
		Table.setColumn("Total Colocolized With Organelle Far", zero_column, results_table_name);
	}
	Table.setColumn("Total Not Colocolized with Organelle", tot_no_colo, results_table_name);
	//Add spot's data (coordinates, intensity and colocalization) per cell
	for (cell = 0; cell < nCells; cell ++) {
		spot_data_per_cell[cell] = "[" + spot_data_per_cell[cell] + "]";
	}
	Table.setColumn("Spots Coordinates Intensity and Colocalization (Far Near or Not Colocolized)", spot_data_per_cell, results_table_name);
	//Add Organelle signal cover per plane columns
	for (column = 0; column < ncross_section; column++) {
		Table.setColumn("Organelle Signal Size Cross Section " + (column + 1), Table.getColumn("Organelle Signal Size Cross Section " + (column + 1), ER_signal_size_name), results_table_name);
	}
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
	close(ER_signal_size_name);
	close(smFISH_list[img]);
//End of img loop
}
print ("All images done. Results saved under " + results_path);
