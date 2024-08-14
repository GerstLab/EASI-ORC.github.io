// Input source folder and create new fodler for converted files
Dialog.create("Important");
Dialog.addMessage("Choose images directory:");
Dialog.show();
bioimg_folder = getDirectory("Choose images directory");
tiff_folder = bioimg_folder + "Tiff Images";
image_format = getString("What is you file format?", "stk");
File.makeDirectory(tiff_folder);
// Setup folders for experiment
channel_num = getNumber("How many channels did you image?", 4);
channels_list = newArray(channel_num);
for (channel_index=0; channel_index < channel_num; channel_index++) {
	channel_name = getString("Name of channel " + channel_index + 1 + " (must be part of file name, case sensative):", "");
	channel_folder = tiff_folder + File.separator + channel_name;
	File.makeDirectory(channel_folder);
	channels_list[channel_index] = channel_name;
}
// Convert bioimg files to tiff
file_num = getFileList(bioimg_folder);
setBatchMode(true);
for (file_index=0; file_index < file_num.length; file_index++) {
	if (endsWith(file_num[file_index], image_format)) {
		image_name = bioimg_folder + file_num[file_index];
		open(image_name);
		for (channel_list_index=0; channel_list_index < channels_list.length; channel_list_index++) {
			if (indexOf(file_num[file_index], channels_list[channel_list_index]) != -1) {
				current_channel_name = channels_list[channel_list_index];
				tiff_path = tiff_folder + File.separator + current_channel_name + File.separator + file_num[file_index];
				saveAs("tiff", tiff_path);
				close();
			}
		}
	}
}