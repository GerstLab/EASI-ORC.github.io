//Getting images directory from user and saving lists of files
Dialog.create("Important");
Dialog.addMessage("Choose proximity marker masks directory:");
Dialog.show();
dapi_folder = getDirectory("Choose proximity marker masks directory:");
Dialog.create("Important");
Dialog.addMessage("Choose Organelle masks directory:");
Dialog.show();
er_folder = getDirectory("Choose Organelle masks directory:");

//Get list of files
dapi_list = getFileList(dapi_folder);
er_list = getFileList(er_folder);

//Create new folders to save sub-segmented masks
save_path = File.getParent(er_folder);
File.makeDirectory(save_path + File.separator + "Near Masks");
File.makeDirectory(save_path + File.separator + "Distant Masks");

dilate_num = parseInt(getString("How many times should the proximity marker mask be dilated (test to see what gets you the best results for your images)?", "3"));
setBatchMode(true);
//Go over images 
for (file = 0; file < dapi_list.length; file++) {
	open(dapi_folder + dapi_list[file]);
	dapi_mask = getTitle();
	// Dilate signal as many times as user chooses
	selectWindow(dapi_mask);
	run("Options...", "iterations="+ dilate_num + " count=1 pad edm=32-bit do=Dilate stack");
	open(er_folder + er_list[file]);
	selectWindow(er_list[file]);
	main_er_mask = getTitle();
	imageCalculator("Subtract create stack", main_er_mask, dapi_mask);
	saveAs("tif", save_path + File.separator + "Distant Masks" + File.separator + main_er_mask + "_distant");
	cER_mask = getTitle();
	imageCalculator("Subtract create stack", main_er_mask, cER_mask);
	saveAs("tif", save_path + File.separator + "Near Masks" + File.separator + main_er_mask + "_near" );
	run("Close All");
}
