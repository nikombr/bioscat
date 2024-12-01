% Prepare a txt that we can input in comsol
protein_structure = 'demoleus2x2'; % Retinin2x2 demoleus2x2

% Load nanostructure
dir = "../../../../../../../work3/s194146/bioscatdata";
xs = load(sprintf("%s/Data/nanostructures/2D/%s_x_1000.txt",dir,protein_structure));
fs = load(sprintf("%s/Data/nanostructures/2D/%s_f_1000.txt",dir,protein_structure));

%figure, plot(xs,fs), grid on, axis equal
surface = [xs' fs'];
save(sprintf('%s/Data/comsol/surface_%s.txt',dir,protein_structure),'surface','-ascii')