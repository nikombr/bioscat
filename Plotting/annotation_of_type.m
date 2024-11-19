function annotation_of_type(scenario,protein_structure,type_anno)

if nargin < 3

    dim = [0.07, 0.03, 0.03, 0.05];
        
    if scenario == 1
        annotation('rectangle',dim,'FaceColor','blue','FaceAlpha',.2,'LineWidth', 1.2)
    elseif scenario == 2
        annotation('rectangle',dim,'FaceColor','red','FaceAlpha',.2,'LineWidth', 1.2)
    end
    annotation('textbox', dim, 'String', sprintf('\\textbf{%d}',scenario), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', 'Interpreter', 'latex','FitBoxToText','on','fontsize',13);
    
    dim = [0.04, 0.03, 0.03, 0.05];
    if strcmp(protein_structure,'demoleus2x2')
        annotation('rectangle',dim,'FaceColor','green','FaceAlpha',.2,'LineWidth', 1.2)
        annotation('textbox', dim, 'String', '\textbf{N}', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', 'Interpreter', 'latex','FitBoxToText','on','fontsize',13);
    elseif strcmp(protein_structure,'Retinin2x2')
        annotation('rectangle',dim,'FaceColor','magenta','FaceAlpha',.2,'LineWidth', 1.2)
        annotation('textbox', dim, 'String', '\textbf{A}', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', 'Interpreter', 'latex','FitBoxToText','on','fontsize',13);
    end

else
    
    dim = [0.04, 0.03, 0.03, 0.05];

    if strcmp(type_anno,'scenario')

        if scenario == 1
            annotation('rectangle',dim,'FaceColor','blue','FaceAlpha',.2,'LineWidth', 1.2)
        elseif scenario == 2
            annotation('rectangle',dim,'FaceColor','red','FaceAlpha',.2,'LineWidth', 1.2)
        end
        annotation('textbox', dim, 'String', sprintf('\\textbf{%d}',scenario), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'EdgeColor', 'none', 'Interpreter', 'latex','FitBoxToText','on','fontsize',13);


    elseif strcmp(type_anno,'protein_type')

            if strcmp(protein_structure,'demoleus2x2')
                annotation('rectangle',dim,'FaceColor','green','FaceAlpha',.2,'LineWidth', 1.2)
                annotation('textbox', dim, 'String', '\textbf{N}', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'EdgeColor', 'none', 'Interpreter', 'latex','FitBoxToText','on','fontsize',13);
            elseif strcmp(protein_structure,'Retinin2x2')
                annotation('rectangle',dim,'FaceColor','magenta','FaceAlpha',.2,'LineWidth', 1.2)
                annotation('textbox', dim, 'String', '\textbf{A}', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'EdgeColor', 'none', 'Interpreter', 'latex','FitBoxToText','on','fontsize',13);
            end

    end

end