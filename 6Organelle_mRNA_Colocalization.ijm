//get the location of cells masks images
Dialog.create("Cells masks location");
Dialog.addMessage("Please input the cells mask images directory:");
Dialog.show();
cell_dir = getDirectory("Cells mask images directory:");
cell_list = getFileList(cell_dir);
//get the location of smFISH csv files
Dialog.create("smFISH csv files location");
Dialog.addMessage("Please input the csv files directory:");
Dialog.show();
smFISH_dir = getDirectory("smFISH csv files directory (created by RSFISH):");
smFISH_list = getFileList(smFISH_dir);
//ask user if ER is segmented and get proper mask images
Dialog.create("ER Sub-segmentation");
Dialog.addChoice("Was ER segmented to nEr and cER using DAPI (Select 'No' if you have one organelle mask)?", newArray("Yes", "No"));
Dialog.show();
ER_seg_check = Dialog.getChoice();
if (ER_seg_check == "Yes") {
	Dialog.create("nER masks location");
	Dialog.addMessage("Please input the perinuclear ER mask images directory:");
	Dialog.show();
	nER_dir = getDirectory("Perinuclear ER mask images directory:");
	nER_list = getFileList(nER_dir);
	Dialog.create("cER masks location");
	Dialog.addMessage("Please input the cortical ER mask images directory:");
	Dialog.show();
	cER_dir = getDirectory("Cortical ER mask images directory:");
	cER_list = getFileList(cER_dir);
	//get the location of DAPI mask images 
	Dialog.create("DAPI masks location");
	Dialog.addMessage("Please input the DAPI mask images directory:");
	Dialog.show();
	DAPI_dir = getDirectory("DAPI mask images directory:");
	DAPI_list = getFileList(DAPI_dir);
}
if (ER_seg_check == "No") {
	Dialog.create("Organelle masks location");
	Dialog.addMessage("Please input the Organelle mask images directory:");
	Dialog.show();
	whole_ER_dir = getDirectory("Organelle mask images directory:");
	whole_ER_list = getFileList(whole_ER_dir);
}
Dialog.create("Results tables location");
Dialog.addMessage("Please input the directory you wish results tables will be saved:");
Dialog.show();
results_path = getDirectory("Results Tables Location:");
smFISH_signal_cutoff = parseFloat(getString("Approximate radius of an smFISH signal (in pixels:", "1.5"));
ER_signal_min = parseFloat(getString("Minimum amount of ER signal within a cell:", "20"));
ER_signal_max = parseFloat(getString("Maximum amount of ER signal within a cell:", "80"));


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
		// Go over non-segmented ER image and create a file from each slice in a new directory
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
	 	if (ER_seg_check == "Yes") {
	 		tot_colo_nER = newArray(nCells);
	 		tot_colo_cER = newArray(nCells);	
	 	}
	 	if (ER_seg_check == "No") {
	 		tot_colo_ER = newArray(nCells);
	 	}
	 	tot_no_colo = newArray(nCells);
	 	ER_signal_size = newArray(nCells);
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
			if (ER_seg_check == "Yes") {
				//Register DAPI signal in each cell (disregard planes where there's no DAPI signal in cell - only for sub-segmented ER)
				open(dapi_slice);
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
					if (Table.get("%Area", roi, cER_signal_table) +  Table.get("%Area", roi, nER_signal_table) < ER_signal_min || Table.get("%Area", roi, cER_signal_table) +  Table.get("%Area", roi, nER_signal_table) > ER_signal_max) {
						ER_signal_size[roi] = 1;
					}
				}
			}
			if (ER_seg_check == "No") {
				//Register % of organelle signal within cell (disregard planes where signal is too great/small)
				selectWindow(current_ER);
				roiManager("select", Array.getSequence(nCells));
				roiManager("measure");
				ER_signal_table = "ER Signal" + plane;
				IJ.renameResults(ER_signal_table);
				selectWindow(current_ER);
				for (roi = 0; roi < nCells; roi++) {
					if (Table.get("%Area", roi, ER_signal_table) < ER_signal_min || Table.get("%Area", roi, ER_signal_table) > ER_signal_max) {
						ER_signal_size[roi] = 1;
					}	
					
				}
			}
			//Go over mRNAs in this slice, find each's cell and measure colocalization
			num_of_mRNA = Table.size(smFISH_list[img]);
			for (mRNA = smFISH_counter; mRNA < num_of_mRNA; mRNA++) {
				//Check "z" column for spot's slice. If z is smaller than current slice + 1, the spot is in current slice.
				if (Table.get("z", mRNA, smFISH_list[img]) >= (plane + 1)) {
					break;
				}
				smFISH_counter += 1;
				//Check for ER, DAPI and edge limits for each cell
				for (cell = 0; cell < nCells; cell++) {
					if (edge_cells_list[cell] == 1 || ER_signal_size[cell] == 1) {
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
					if (Roi.contains(current_X, current_Y)) {
						//Use ability to draw selection circles on organelle mask image - sub-segmented ER
						if (ER_seg_check == "Yes") {
							selectWindow(current_cER);
							//Create circle selection with radius of cutoff given by user, around smFISH spot pixel
							makeOval(current_X - smFISH_signal_cutoff, current_Y - smFISH_signal_cutoff, smFISH_signal_cutoff * 2, smFISH_signal_cutoff * 2);
							run("Clear Results");
							run("Set Measurements...", "mean redirect=None decimal=2");
							run("Measure");
							selectWindow(current_nER);
							//Create circle selection with radius of cutoff given by user, around smFISH spot pixel
							makeOval(current_X - smFISH_signal_cutoff, current_Y - smFISH_signal_cutoff, smFISH_signal_cutoff * 2, smFISH_signal_cutoff * 2);
							run("Measure");
							//Check if there's ER signal near the smFISH spot
							if (getResult("Mean", 0) != 0 && getResult("Mean", 1) != 0) {
								if (getResult("Mean", 0) > getResult("Mean", 1)) {
									tot_colo_cER[cell]++;
									continue;											
								}
								else if (getResult("Mean", 0) <= getResult("Mean", 1)) {
									tot_colo_nER[cell]++;
									continue;
								}
							}
							else if (getResult("Mean", 0) != 0) {
							tot_colo_cER[cell]++;
							continue;										
							}
							else if (getResult("Mean", 1) != 0) {
								tot_colo_nER[cell]++;
								continue;
							}
							else {
								tot_no_colo[cell]++;
								continue;
							}
						}
						//Use ability to draw selection circles on organelle mask image - single mask organelle
						if (ER_seg_check == "No") {
							selectWindow(current_ER);
							//Create circle selection with radius of cutoff given by user, around smFISH spot pixel
							makeOval(current_X - smFISH_signal_cutoff, current_Y - smFISH_signal_cutoff, smFISH_signal_cutoff * 2, smFISH_signal_cutoff * 2);
							run("Clear Results");
							run("Set Measurements...", "mean redirect=None decimal=2");
							run("Measure");
							//Check if there's ER signal near the smFISH spot
							if (getResult("Mean", 0) != 0) {
								tot_colo_ER[cell]++;
								continue;
							}
							else {
								tot_no_colo[cell]++;
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
	}
	//Sum and save data for img
	results_table_name = cell_list[img] + "Colocalization.csv";
	Table.create(results_table_name);
	Table.setColumn("Cell #", Array.getSequence(tot_mRNA_per_cell.length), results_table_name);
	Table.setColumn("Total mRNA per Cell", tot_mRNA_per_cell, results_table_name);
	if (ER_seg_check == "Yes") {
		Table.setColumn("Total Colocolized With nER", tot_colo_nER, results_table_name);
		Table.setColumn("Total Colocolized With cER", tot_colo_cER, results_table_name);	
	}
	if (ER_seg_check == "No") {
		Table.setColumn("Total Colocolized With Organelle", tot_colo_ER, results_table_name);
	}
	Table.setColumn("Total Not Colocolized with Organelle", tot_no_colo, results_table_name);
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
	close(smFISH_list[img]);
//End of img loop
}
print ("All images done. Results saved under " + results_path);
