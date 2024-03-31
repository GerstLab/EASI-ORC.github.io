Dialog.create("Important");
Dialog.addMessage("Please make sure to use RS-FISH manually for each probe you're using, to identify proper threshold.");
Dialog.show();
//Ask user for cells mask images location
Dialog.create("Cells masks location");
Dialog.addMessage("Please input the cells mask images directory:");
Dialog.show();
cell_dir = getDirectory("Cells mask images directory:");
cell_list = getFileList(cell_dir);
//Ask user for location of FISH images and threshold for RSFISH
Dialog.create("smFISH images location");
Dialog.addMessage("Please input the smFISH (probe) images directory:");
Dialog.show();
smFISH_dir = getDirectory("smFISH images directory (probe channel):");
RSFISH_Threshold = parseFloat(getString("Threshold for RS FISH:", "0.00180"));
smFISH_list = getFileList(smFISH_dir);
File.makeDirectory(smFISH_dir + "smFISH Spots");
//Go over each image in the folder and create ROI files with all the spots identified
for (img = 0; img < smFISH_list.length; img++) {
	FISH_img_path = smFISH_dir + smFISH_list[img];
	cells_mask_path = cell_dir + cell_list[img];
	if (!File.isDirectory(FISH_img_path)) {
		//Use cell masks to remove all areas outside of cells in smFISH image
		open(cells_mask_path);
		//Dilate to improve edge dot identification and use Watershed to seperate close cells
		setOption("BlackBackground", false);
		run("Dilate", "stack");
		run("Watershed", "stack");
		//Analyze all cells and add to roi manager
		run("Analyze Particles...", "exclude add slice");
		run("Close All");
		open(FISH_img_path);
		title = getTitle();
		//Calculate mean intensity of image, to use as threshold
		roiManager("deselect");
		close("ROI Manager");
		selectWindow(title);
		run("Select All");
		num_of_pix = getWidth() * getHeight();
		sum_of_intesities = 0;
		run("Set Measurements...", "integrated redirect=None decimal=5");
		for (slice = 1; slice <= nSlices; slice++) {
			setSlice(slice);
			run("Measure");
			sum_of_intesities += getResult("RawIntDen", slice - 1);
		}
		close("Results");
		//Run RSFISH
		run("RS-FISH", "image=[" + title + "] mode=Advanced anisotropy=1.0000 robust_fitting=RANSAC compute_min/max use_anisotropy spot_intensity=[Linear Interpolation] sigma=1.50000 threshold=" + RSFISH_Threshold + " support=3 min_inlier_ratio=0.10 max_error=1.50 spot_intensity_threshold=0 background=[No background subtraction] background_subtraction_max_error=0.05 background_subtraction_min_inlier_ratio=0.10 results_file=[] num_threads=16 block_size_x=128 block_size_y=128 block_size_z=16");
		//Save results table
		selectWindow("smFISH localizations");
		table_name = title + " smFISH localizations.csv";
		table_path = smFISH_dir + "smFISH Spots" + File.separator + title + " smFISH localizations.csv";
		if (Table.size >= 1) {
			Table.sort("z");
			saveAs("Results", table_path);
		}
		else {
			//setResult("x", nResults, "No spots");
			saveAs("Results", table_path);
		}
		//Close all windows before starting next image
		close(table_name);
		run("Close All");
		close("Log");
	}
}
