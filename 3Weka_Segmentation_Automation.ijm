//Getting image directory from user and creating probability maps folder
bioimg_folder = getDirectory("Input image directory:");
prob_map_folder = bioimg_folder + File.separator + "Probabiliy Maps";
File.makeDirectory(prob_map_folder);
//Getting trained calassifier file location from user
classifier_path =  File.openDialog("Input classifier location:");

//Starting time
startTime = getTime();

//Get list of images
file_list = getFileList(bioimg_folder);
setBatchMode(true);
for (file = 0; file < file_list.length; file++) {
	//Make sure file is not a folder (won't go into fsubfolders)
	image_path = bioimg_folder + file_list[file];
	//Create folders for seperated slices and maps
	seperated_folder = bioimg_folder + "seperated_planes";
	File.makeDirectory(seperated_folder);
	temp_maps_folder = bioimg_folder + "temp_maps";
	File.makeDirectory(temp_maps_folder);
	if (File.isFile(image_path)) {
	 	//Go over files and create a file from each slice in a new directory
	 	open(image_path);
	 	run("Image Sequence... ", "dir=[" + seperated_folder + "] format=TIFF digits=2 use");
		run("Close All");
		//Go over new folder, and use classifier to output probability map
		plane_list = getFileList(seperated_folder);
		for (plane = 0; plane < plane_list.length; plane++) {
			setBatchMode(false);
			open(seperated_folder + File.separator + plane_list[plane]);
			run("Trainable Weka Segmentation");
			while (!isOpen("Trainable Weka Segmentation v3.3.2")) {
				wait(100);
			}
			selectWindow("Trainable Weka Segmentation v3.3.2");
			call("trainableSegmentation.Weka_Segmentation.loadClassifier", classifier_path);
			call("trainableSegmentation.Weka_Segmentation.getProbability");
			while (!isOpen("Probability maps")) {
				wait(100);
			}
			selectWindow("Probability maps");
			run("Delete Slice");
			saveAs("tiff", temp_maps_folder + File.separator + plane_list[plane]);
			//Force garbage collection (important for large images)
			run("Close All");
           	call("java.lang.System.gc")
		}
		setBatchMode(true);
		//Create a multi-plane image from probability maps of all planes
		map_list = getFileList(temp_maps_folder);
		for (map = 0; map < map_list.length; map++) {
			map_path = temp_maps_folder + File.separator + map_list[map];
			open(map_path);
			}
		// Save maps as stack
		run("Images to Stack", "name=" + file_list[file] + "");
		saveAs("tiff", prob_map_folder + File.separator + file_list[file] + "_map");
		run("Close All");
		//Delete uneeded files and folders
		sepereated_list = getFileList(seperated_folder);
		for (sep = 0; sep < sepereated_list.length; sep++) {
			img_to_del = seperated_folder + File.separator + sepereated_list[sep];
			File.delete(img_to_del);
		}
		maps_list = getFileList(temp_maps_folder);
		for (dmap = 0; dmap < maps_list.length; dmap++) {
			map_to_del = temp_maps_folder + File.separator + maps_list[dmap];
			File.delete(map_to_del);
		}
		File.delete(seperated_folder + File.separator);
		File.delete(temp_maps_folder + File.separator);
	}
}
//Print elapsed time
estimatedTime = (getTime() - startTime) * 0.001;
IJ.log( "** Finished processing folder in " + estimatedTime + " s **" );