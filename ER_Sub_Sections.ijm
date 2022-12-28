//Getting images directory from user and saving lists of files
dapi_folder = getDirectory("Choose DAPI probability maps directory:");
er_folder = getDirectory("Choose ER probability maps directory:");

dapi_list = getFileList(dapi_folder);
er_list = getFileList(er_folder);

dilate_num = parseInt(getString("How many times should the DAPI mask file be dilated (test to see what gets you the best results for your images)?", "3"));
setBatchMode(true);
//Go over images 
for (file = 0; file < dapi_list.length; file++) {
	open(dapi_folder + dapi_list[file]);
	//Convert to mask file
	run("Make Binary", "method=Minimum background=Light calculate create");
	dapi_mask = getTitle();
	//Save DAPI mask files
	File.makeDirectory(dapi_folder + "DAPI Masks");
	selectWindow(dapi_mask);
	saveAs("tif", dapi_folder + "DAPI Masks" + File.separator + dapi_mask);
	// Dilate signal as many times as user chooses
	selectWindow(dapi_mask);
	run("Options...", "iterations="+ dilate_num + " count=1 pad edm=32-bit do=Dilate stack");
	open(er_folder + er_list[file]);
	selectWindow(er_list[file]);
	run("Make Binary", "method=Li background=Light calculate create");
	main_er_mask = getTitle();
	imageCalculator("Subtract create stack", main_er_mask, dapi_mask);
	File.makeDirectory(er_folder + "cER Masks");
	saveAs("tif", er_folder + "cER Masks" + File.separator + main_er_mask + "_cER");
	cER_mask = getTitle();
	imageCalculator("Subtract create stack", main_er_mask, cER_mask);
	File.makeDirectory(er_folder + "nER Masks");
	saveAs("tif", er_folder + "nER Masks" + File.separator + main_er_mask + "_nER" );
	run("Close All");
}
