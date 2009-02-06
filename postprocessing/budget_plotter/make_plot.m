function zzz = make_plot (nnn,var1,grid1,marker1,title1,units1,maxdim,zbottom,ztop)

evalstring = ['top',nnn,' = subplot (2,5,',nnn,');'];
eval(evalstring);
%setstring = ['top',nnn,',''FontSize'',8'];
%setstring
%eval(setstring)
set(gca,'FontSize',8);

hold on

if (or((nnn == '1'),(nnn == '6')))
    plot(var1,grid1,marker1,'Linewidth',2);
elseif (or((nnn == '4'),(nnn == '9')))
    plot(var1,grid1,marker1,'Color',[1,0.5,0]);
else
    plot(var1,grid1,marker1);
end

line([0 0],[zbottom ztop],'Color','k','Linestyle','-');

title([title1]);

xlabel(units1);
xlim([maxdim*-1 maxdim]);
ylim([zbottom ztop]);
hold off