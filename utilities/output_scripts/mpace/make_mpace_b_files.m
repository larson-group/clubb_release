function [] = make_mpace_b_files();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make_mpace_b_files()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Michael Falk, September 2006-March 2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters: none
% Input files: mpace_b GrADS files (.ctl and .dat)
% Output parameters: none
% Output files:  
%
% Requires: two Michael Falk utility scripts:
%           header_read.m
%           read_grads_hoc_endian.m
%
% Reference: http://science.arm.gov/workinggroup/cpm/scm/scmic5/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program makes the required output for submission to ARM for the
% mpace_b (Mixed-Phase Arctic Cloud Experiment) case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set up input files and timestep
scm_path = ['/home/mjfalk/hoc/hoc_mpace_b_0047k/'];
smfile   = 'mpace_b_zt.ctl';

t = 0:30:720;
t(1) = 1;
sizet = size(t);
sizet = max(sizet);

% Load mpace_b data from mass grid files
[filename,nz,z,ntimesteps,numvars,list_vars] = header_read([scm_path,smfile]);

for i=1:numvars
    for timestep = 1:sizet-1
        stringtoeval = [list_vars(i,:), ' = read_grads_hoc_endian([scm_path,filename],''ieee-le'',nz,t(timestep),t(timestep+1),i,numvars);'];
        eval(stringtoeval);
        str = list_vars(i,:);
            sprintf('%s','Pressure File');
            sprintf('%d %d',nz,sizet-1');
        for j=1:nz
            arraydata(j,timestep) = eval([str,'(j)']);
        end
    end
    eval([strtrim(str),'_array = arraydata;']);
end

% Set up constants
g0 = 9.8;
p0 = 1e5;
R  = 287.04;
Cp = 1004.67;
Lv = 2.5e6;
Ni = 2000; % ice concentration per liter for simplified ice scheme

% Set derived variables
exner_array = (p_array./p0).^(R/Cp);
T_array = (thlm_array .* exner_array) + (Lv/Cp).*rcm_array;
rvm_array = rtm_array - rcm_array;
RH_array = rvm_array ./ rsm_array;

rsnowm_array = (Ni * m_array) ./ (rhot_array);
snc_array = (thlm_array .* 0) + Ni;

size_snow_array = size(rsnowm_array);
snow_width = size_snow_array(1);
snow_length = size_snow_array(2);
for i = 1:snow_width
    for j = 1:snow_length
        if (rsnowm_array(i,j) > .000001)
            hydro_array(i,j) = 1;
        else
            hydro_array(i,j) = 0;
        end
        
        if (cf_array(i,j) == 1)
            hydro_array(i,j) = 1;
        end
    end
end


% Set up list of variables and output file parameters
variablelist = ['p          ';'T          ';'rvm        ';'RH         ';
                'rcm        ';'rsnowm     ';'cf         ';'um         ';
                'vm         ';'radht      ';'radht_SW   ';'radht_LW   ';
                'hydro      ';'snc        ';];

filenumber   = [0;2;3;4;5;8;10;11;12;19;20;21;31;36;];
                
descriptions = ['Pressure (mb)                  ';
                'Temperature (K)                ';
                'Water vapor mixing ratio (g/kg)';
                'Relative humidity ()           ';
                'Cloud water mixing ratio (g/kg)';
                'Snow mixing ratio (g/kg)       ';
                'Cloud fraction ()              ';
                'U-Wind (m/s)                   ';
                'V-Wind (m/s)                   ';
                'Radiative Heating Rate (K/day) ';
                'Shortwave Heating Rate (K/day) ';
                'Longwave Heating Rate (K/day)  ';
                'Hydrometeor Fraction ()        ';
                'Snow number concentration (#/l)'; ];

% This scales each variable from how it appears in GrADS output to how
% ARM wants it submitted.
multfactor   = [.01; 1; 1000; 1; 1000; 1000; 1; 1; 1; 86400; 86400; 86400; 1; .001;];

formatting   = ['%07.1f';'%07.2f';'%07.3f';'%06.3f';'%07.4f';'%07.4f';'%06.3f';
                '%07.2f';'%07.2f';'%07.2f';'%07.2f';'%07.2f';'%06.3f';'%07.1f'];

% Each variable is written to a separate ASCII text file.
% Half-hourly profiles:
sizevars = size(variablelist);
for v=1:sizevars(1)
    % open file
    str = strtrim(variablelist(v,:))
    string = ['fid=fopen(''b1_p',num2str(filenumber(v)),'_UWM'',''w'')'];
    eval(string);
   
    % write header and max/min
    string = ['#',str,'  UWM  ',descriptions(v,:),'  mjfalk  19 Oct 2006'];
    fprintf(fid,'%s\n',string);
    fprintf(fid,'%d %d\n',nz,sizet-1);
    
    string   = ['max(max(',str,'_array))'];
    maxvalue = eval(string);
    string   = ['min(min(',str,'_array))'];
    minvalue = eval(string);
    string = ['fprintf(fid,''',formatting(v,:),'  ',formatting(v,:),'\n'',',num2str(maxvalue*multfactor(v)),',',num2str(minvalue*multfactor(v)),')';]
    eval(string);
    
    % write vertical column of pressures
    for j=1:nz
        fprintf(fid,'%07.1f %s',p(j,1)/100,' '); % Converting from Pa to mb
        if (mod(j,10) == 0)
            fprintf(fid,'%s\n','');
        end
    end

    fprintf(fid,'%s\n','');

    % write list of output times
    for timestep=1:sizet-1
        fprintf(fid,'%7.2f %s',(1/60)*((t(timestep+1)-t(timestep))/2 + t(timestep)),' ');
        if (mod(timestep,10) == 0)
            fprintf(fid,'%s\n','');
        end
    end

    fprintf(fid,'%s\n','');

    % write data in half-hour averages
    for j=1:nz
        for timestep=1:sizet-1
            string = ['value = ',str,'_array(',num2str(j),',',num2str(timestep),');'];
            eval(string);
            alldata(j,timestep) = value;
            if (strcmp(str,'snc'))
                string = ['fprintf(fid,''',formatting(v,:),' %s'',value*multfactor(v)*hydro_array(',num2str(j),',',num2str(timestep),'),'' '');'];
            else
                string = ['fprintf(fid,''',formatting(v,:),' %s'',value*multfactor(v),'' '');'];
            end
             eval(string);
            if (mod(timestep,10) == 0)
                fprintf(fid,'%s\n','');
            end
        end
        fprintf(fid,'%s\n','');
    end

    % Plot everything and step through with keyboard command.
    if (v > 1)
        for plotloop=1:(sizet-1)
            clf
            plot(alldata(:,plotloop).*multfactor(v),p(:,1))
            set(gca,'YDir','reverse')
            title([descriptions(v,:),' at time ',num2str(plotloop)])
            keyboard
        end
    end
end


% Make files full of zeroes, for the variables we're not submitting.
variablelist = ['rim        ';'rrm        ';'rgm        ';'q1c        ';
                'q2         ';'r_eff_l    ';'r_eff_i    ';'ncm        ';
                'nim        ';'gnc        '];

filenumber   = [6;7;9;13;14;32;33;34;35;37;];
                
descriptions = ['Cloud ice mixing ratio (g/kg)             ';
                'Rain water mixing ratio (g/kg)            ';
                'Graupel mixing ratio (g/kg)               ';
                'Apparent heat source (K/day)              ';
                'Apparent moisture sink (K/day)            ';
                'Liquid cloud droplet effective radius (um)';
                'Ice crystal effective size (um)           ';
                'Cloud droplet number (#/cm^3)             ';
                'Cloud ice number concentration (#/l)      ';
                'Graupel number concentration (#/l)        '; ];

multfactor   = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1;];

formatting   = ['%07.4f';'%07.4f';'%07.4f';'%07.2f';'%07.2f';
                '%07.2f';'%07.2f';'%07.1f';'%07.1f';'%07.1f'];

sizevars = size(variablelist);
for v=1:sizevars(1)
    str = strtrim(variablelist(v,:))
    string = ['fid=fopen(''b1_p',num2str(filenumber(v)),'_UWM.not_submitted'',''w'')'];
    eval(string);
end




% Load mpace_b data from momentum grid files
swfile   = 'mpace_b_zm.ctl';
[filename2,nz2,z2,ntimesteps2,numvars2,list_vars2] = header_read([scm_path,swfile]);

for i=1:numvars2
    for timestep = 1:sizet-1
        stringtoeval = [list_vars2(i,:), ' = read_grads_hoc_endian([scm_path,filename2],''ieee-le'',nz2,t(timestep),t(timestep+1),i,numvars2);'];
        eval(stringtoeval);
        str = list_vars2(i,:);
            sprintf('%s','Pressure File');
            sprintf('%d %d',nz2,sizet-1');
        for j=1:nz2
            arraydata(j,timestep) = eval([str,'(j)']);
        end
    end
    eval([strtrim(str),'_array = arraydata;']);
end

% Create three time series files:
for v=1:3
    idstr = ['Group ',num2str(v)];
    string = ['fid=fopen(''b1_t',num2str(v),'_UWM'',''w'');'];
    eval(string);
   
    string = ['#',num2str(idstr),'  UWM  mjfalk  19 Oct 2006'];
    fprintf(fid,'%s\n',string);
    fprintf(fid,'%d\n',nz);
    
    for timestep=1:sizet-1
        listoftimes(timestep) = (1/60)*((t(timestep+1)-t(timestep))/2 + t(timestep));
    end

    % Set up the variables for file 1
    if (v==1)
        field1 = listoftimes;
        field2 = (T_array(2,:) * zeros + ones) * -888.89; % SST (K) F7.2
        field3 = Cp*T_array(2,:) + g0*z(2); % Dry static energy (kJ/kg) (factor .001) F7.2
        field4 = rvm_array(2,:); %Water vapor mixing ratio (g/kg) (factor 1000) F6.2
        field5 = field3 + Lv.*field4; % Moist static energy (kJ/kg) (factor .001) F7.2
        field6 = um_array(2,:); % u wind (m/s) F7.2
        field7 = vm_array(2,:); % v wind (m/s) F7.2
        field8 = rhom_array(1,:) .* Cp .* wpthlp_array(1,:); % sfc turb flux of sens. heat (W/m^2) F6.1
        field9 = rhom_array(1,:) .* Lv .* wprtp_array(1,:); %sfc turb flux of latent heat (W/m^2) F6.1
        field10= rhom_array(1,:) .* upwp_array(1,:); %sfc turb flux of momentum (N/m^2) F8.4
        field11= rhom_array(1,:) .* vpwp_array(1,:); %sfc turb flux of momentum (N/m^2) F8.4

        multfactor   = [1; 1; .001; 1000; .001; 1; 1;
                        1; 1; 1; 1;];

        formatting   = ['%07.2f';'%07.2f';'%07.2f';'%06.2f';'%07.2f';
                        '%07.2f';'%07.2f';'%06.1f';'%06.1f';'%08.4f';'%08.4f';];


        % plot everything and step through with keyboard command
        for lup = 1:11
            eval(['plot(30:30:720,field',num2str(lup),'*multfactor(',num2str(lup),'))'])
            keyboard
        end
    
    end

    % Set up the variables for file 2
    if (v==2)
        field1 = listoftimes; % F7.1, has to be
        field2 = (T_array(2,:) * zeros + ones) * -888.888888; % Sfc downwelling SW (W/m^2) F7.1
        field3 = field2; % Sfc upwelling SW (W/m^2) F7.1
        field4 = field2; % Sfc downwelling IR (W/m^2) F6.2
        field5 = field2; % Sfc upwelling IR (W/m^2) F6.2
        field6 = field2; % u wind (m/s) F7.2
        field7 = field2; % v wind (m/s) F7.2
        field8 = field2; % sfc turb flux of sens. heat (W/m^2) F6.1
        field9 = cf_array(1,:) .* zeros + ones; %sfc turb flux of latent heat (W/m^2) F6.1

        numberoftimes = size(listoftimes);
        numberoftimes = max(numberoftimes);

        precwat(numberoftimes) = 0;
        lwp(numberoftimes) = 0;
        cip(numberoftimes) = 0;  %% THIS ACTUALLY IS VERTICALLY INTEGRATED SNOW.
        
        for zz=1:nz-1
            for tt=1:numberoftimes
                rho = rhom_array(zz,tt);
                q   = 0.5*(rvm_array(zz,tt)+rvm_array(zz+1,tt));
                qc  = 0.5*(rcm_array(zz,tt)+rcm_array(zz+1,tt));
                qi   = 0.5*(rsnowm_array(zz,tt)+rsnowm_array(zz+1,tt));
                dz  = z(zz+1)-z(zz);
                precwat(tt) = precwat(tt) + (rho * q * dz);
                lwp(tt)     = lwp(tt) + (rho * qc * dz);
                cip(tt)     = cip(tt) + (rho * qi * dz);
                snowrate(tt) = 0.5*(rsnowm_array(1,tt));
            end
        end

        field10= precwat(:)'; %sfc turb flux of momentum (N/m^2) F8.4
        field11= lwp(:)'; %sfc turb flux of momentum (N/m^2) F8.4
        field12 = field2 * zeros;

        multfactor   = [1; 1; 1; 1;
                        1; 1; 1; 1;
                        1; 1; 1; 1;];

        formatting   = ['%07.2f';'%07.2f';'%07.2f';'%06.2f';'%07.2f';
                        '%07.2f';'%07.2f';'%06.1f';'%06.1f';'%08.4f';'%08.4f';'%08.4f';];
        % plot everything and step through with keyboard command
        for lup = 1:12
            eval(['plot(30:30:720,field',num2str(lup),'*multfactor(',num2str(lup),'))'])
            keyboard
        end
    end

    
    % Set up the variables for file 3
    if (v==3)
        field1 = listoftimes;
        field2 = T_array(2,:) * zeros; % Vertically integrated rain
        field3 = cip(:)'; % Dry static energy (kJ/kg) (factor .001) F7.2
        field4 = field2; % Vertically integrated graupel
        field5 = field2; % Surface rain rate

        u_T_cm_sfc = max(u_T_cm_array)
        field6= rsnowm_array(2,:) .* rhot_array(2,:) .* (u_T_cm_sfc(:)'./100); %sfc turb flux of momentum (N/m^2) F8.4

        multfactor(:) = 0;
        multfactor   = [1; 1; 1; 1; 1; 86400;];
        formatting   = [' %07.2f';'%010.3f';'%010.3f';'%010.3f';' %07.2f';' %07.2f';];

        % plot everything and step through with keyboard command
        for lup = 1:6
            eval(['plot(30:30:720,field',num2str(lup),'*multfactor(',num2str(lup),'))'])
        end
    end
    fprintf(fid,'%s\n','');

    % Find maxima and minima
    fprintf(fid,'%s ','+999.9');
    for loop=1:size(multfactor)
        string = ['fprintf(fid,''',strtrim(formatting(loop,:)),' '',',num2str(multfactor(loop)),'*max(field',num2str(loop),'));']
        eval(string)
    end
    fprintf(fid,'%s\n','');

    fprintf(fid,'%s ','-999.9');
    for loop=1:size(multfactor)
        string = ['fprintf(fid,''',strtrim(formatting(loop,:)),' '',',num2str(multfactor(loop)),'*min(field',num2str(loop),'));']
        eval(string)
    end
    fprintf(fid,'%s\n','');

    
    % Print the fields to the files
    for timestep=1:sizet-1
        for nfield = 1:size(multfactor)
            string = ['value = field',num2str(nfield),'(:,',num2str(timestep),');'];
            eval(string);
            string = ['fprintf(fid,''',strtrim(formatting(nfield,:)),' '',value*multfactor(nfield));'];
            eval(string);
        end
        fprintf(fid,'%s\n','');
    end

end

fclose('all');
