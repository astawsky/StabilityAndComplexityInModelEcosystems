function rw = ggg5(t,y)

load('forformula.mat','extinction_point','number_of_species','k_extinction','r','extinct');
load('forformula.mat','A');
%dpop=@(t,y) r(l).*(1-dot(A(l,:),y)).*y(l); l, number_of_species-l, l
dpop_str="r(%d).*(1-dot(A(%d,:),y)).*y(%d)";
%backup=@(l,r,k_extinction,y) (-r(l))*k_extinction*y(l);
backup_str="-k_extinction*y(%d)";
for l=1:number_of_species
    if y(l)>extinction_point && extinct(l)==0
        if l==1
            some="[";
            some=strcat(some,sprintf(dpop_str,l,l,l),";");
        elseif l==number_of_species
            some=strcat(some," ",sprintf(dpop_str,l,l,l));
            some=strcat(some,"];");
        else
            some=strcat(some," ",sprintf(dpop_str,l,l,l),";");
        end
    else
        extinct(l)=1;
        if l==1
            some="[";
            some=strcat(some,sprintf(backup_str,l,l),";");
        elseif l==number_of_species
            some=strcat(some," ",sprintf(backup_str,l,l));
            some=strcat(some,"];");
        else
            some=strcat(some," ",sprintf(backup_str,l,l),";");
        end
    end
end
save('forformula.mat','extinction_point','number_of_species','k_extinction','r','extinct','A');

 ff=eval(strcat('@(t,y) ',some));
 
 rw=ff(t,y);

end