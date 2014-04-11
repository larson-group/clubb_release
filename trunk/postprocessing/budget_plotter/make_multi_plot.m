function zzz = make_multi_plot (nnn,multi_var,multi_grid,multi_marker,multi_legend,legendpos,title1,units1,maxdim,zbottom,ztop)

evalstring = ['top',nnn,' = subplot (2,5,',nnn,');'];
eval(evalstring);
set(gca,'FontSize',8);

sizearray = size(multi_var);
sizearray = sizearray(1);

string    = [];

for i=1:sizearray-1
    j=num2str(i);
    string = [string,'multi_var(',j,',:),multi_grid(',j,',:),multi_marker(',j,',:),'];
end
i = i + 1;
j=num2str(i);
string = [string,'multi_var(',j,',:),multi_grid(',j,',:),multi_marker(',j,',:)'];

%hold all
hold on
eval(['plot(',string,');']);
line([0 0],[zbottom ztop],'Color','k','Linestyle','-');

title([title1]);

xlabel(units1);
xlim([maxdim*-1 maxdim]);
ylim([zbottom ztop]);
%legend(multi_legend,'Location',legendpos);
legend(multi_legend);
legend('boxoff');

hold off