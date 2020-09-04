
clear all;
close all;

% MFILE_ROOT='../../../mfiles';
% path([MFILE_ROOT '/fun-spool'],path);
% path([MFILE_ROOT '/saclab'],path);

% -- parameter --
%mx=49; my=45; mz=26;
%mx=47; my=90; mz=17;
mx=160;my=106;mz=25;
%fnm_blk='../block.49x45x26.1x1x1.2x2x1.nc'
%fnm_blk='../block.47x90x17.1x1x1.1x1x1.nc'
fnm_blk='../block.160x106x25.1x1x1.1x1x1.nc'

%smot_list=[64];
% smot_list=[2 4 8 16 24];
% damp_list=[2 4 8 16 24];
smot_list=[2];
damp_list=[2];
%zval_list=[ 7 15 25 36 52 76];
%zindx_list=[25 23 21 19 16 12];
zindx_list=[24 22 21 20 16 12];
figlabel={'(a) ','(b) ','(c) ','(d) ','(e) ','(f) '};
% zval_list=[52 ];
% zindx_list=[16];
% res_list={'./result.1th.8x8x24'};
res_list={'./result.1th.4x4x24','./result.1th.6x6x24','./result.1th.8x8x24',...
    './result.1th.10x10x24','./result.1th.12x12x24','./result.1th.14x14x24','./result.1th.16x16x24'};
hscale_list=[26, 40, 53, 66, 79, 92, 106];
%res_list={'./result.2th'};
%res_list={'../Freq123.ENZ.dt.lt.10.model.wg1/result.1th','./result.2th'};
cmp_list={'Vs'};

savefig_flag=1;

% ----------------------------------
num_z=numel(zindx_list);
num_smot=numel(smot_list);
num_damp=numel(damp_list);
num_res=numel(res_list);

X=nc_varget(fnm_blk,'x');
Y=nc_varget(fnm_blk,'y');
Z=nc_varget(fnm_blk,'z');

x=90-reshape(X,[mx,my,mz])/pi*180;
y=reshape(Y,[mx,my,mz])/pi*180;
R=reshape(Z,[mx,my,mz]);
z=(6371-reshape(Z,[mx,my,mz])/1e3);

X=reshape(X,[mx,my,mz]);
Y=reshape(Y,[mx,my,mz]);
% convert to Cartesian coord
[yc,xc,zc]=sph2cart(Y,pi/2-X,R);
xc=xc/1e3; yc=yc/1e3; zc=zc/1e3; %in km

minlon=283.8; maxlon=293.2; minlat=38.2;maxlat=47.8;
%minlon=284; maxlon=292; minlat=39;maxlat=47;

% ----------------------------------
%for nres=1:num_res
%figure('Position',[400 400 800 500]); box on;
for nres=1:num_res
    pnm_result=res_list{nres};
    pnm_fig=[pnm_result '.fig'];
    if ~ isdir(pnm_fig)
       mkdir(pnm_fig);
    end
    %
    %temp
    
    %subplot(2,3,nres);
    %end temp
    %
    for nsmot=1:num_smot
        KSMOTNM=num2str(smot_list(nsmot));
        for ndamp=1:num_damp
            KDAMPNM=num2str(damp_list(ndamp));

            fnm_try=[pnm_result '/try.damp' KDAMPNM '.smot' KSMOTNM '.st0.ev0.lo0.dat'];

            M=load(fnm_try);
            W=reshape(M,mx,my,mz,2);
            %    W=permute(W,[2,1,3,4]);

            W=W*100; %percent

            % z slice
            for ncmp=1:length(cmp_list)
                KCMPNM=cmp_list{ncmp};
                if strcmp(KCMPNM, 'Vp')
                    Vel=squeeze(W(:,:,:,1));
                elseif strcmp(KCMPNM, 'Vs')
                    Vel=squeeze(W(:,:,:,2));
                else
                    error('***ERROR: Wrong component tag, only: Vp and/or Vs are valid!');
                end

                fnm_fig=[pnm_result,'_',cmp_list{ncmp},'_s' num2str(smot_list(nsmot)) 'd' num2str(damp_list(ndamp))];
                figure('Position',[400 400 800 500]);
                %{
                %     temp
                % subplot(2,3,1)
                % pcolor(YS,XS,V);
                % 
                %         shading interp
                %         colormap('jetwr');
                %         %colormap('jet');
                % 
                %         daspect([1,cosd((minlat+maxlat)/2),1]);
                %         xlabel('Longitude, deg.'); ylabel('Latitude, deg.');
                %         axis([minlon maxlon minlat maxlat]);
                %         colorbar
                % 
                %         title([figlabel{1} 'input model'])
                %         caxis([-10,10]);

                % end temp
                %}    


                for nk=1:num_z
                    k=zindx_list(nk);
                    %d=zval_list(nk);
                    V=squeeze(Vel(:,:,k));
                    % remove the mean, 04/26/2011, Y.Shen
                    V=V-mean(mean(V));
                    %
                    XS=squeeze(x(:,:,k));
                    YS=squeeze(y(:,:,k));
                    zs=squeeze(zc(:,:,k));
                    %KDNM=num2str(zval_list(nk));
                    KINM=num2str(zindx_list(nk));
                    Vmax=max(max(abs(V)));

                    fnm_outfile=[cmp_list{ncmp},'_' int2str(nk) '.txt'];
                    Vmax=10;

                    subplot(2,3,nk);
                    %set(hid,'renderer','zbuffer');
                    set(gca,'TickDir','out');

                    %       imagesc(x1d,y1d,V);hold on
                    %	imagesc(y1d,x1d,V);hold on
                    hold on;
                    h=image(YS(1,:),XS(:,1),V);
                    h.CDataMapping='scaled';
                    %pcolor(YS(1,:),XS(:,1),V);

                    %shading flat;
                    colormap('jetwr');
                    colorbar;
                    set(gca,'CLim',[-Vmax,Vmax]);
                    %colormap('jet');
                    
                    state=[];
                    load us_states;
                    for sb=1:length(state)
                        plot(state(sb).polygon(:,1)+360, state(sb).polygon(:,2),'color','k','LineWidth',1);
                    end
                    
                    % xlim(scl_xlim); ylim(scl_ylim);
                    daspect([1,cosd((minlat+maxlat)/2),1]);
                    %xlabel('Longitude, deg.'); ylabel('Latitude, deg.');
                    axis([minlon maxlon minlat maxlat]);
                    text(maxlon-2.5,minlat+1,strcat(num2str(round(z(round(mx/2),round(my/2),k))),' km'));
                    

%                     title(strcat(figlabel{nk}, 'scale: ',num2str(hscale_list(nres)), ...
%                         ' km, z: ',num2str(round(z(round(mx/2),round(my/2),k))),' km'))
                    title(strcat(figlabel{nk}, 'cell size: ',num2str(hscale_list(nres)), ...
                        ' km'));
                    
                    %caxis([-Vmax,Vmax]);
                    box on;
                    hold off;
                end
                if savefig_flag
                    saveas(gca,['./' pnm_fig '/' fnm_fig '.ps'],'psc');
                    pause;
                    close;
                end
            end
        end
    end
    drawnow
end

%plot([55 55],[23,43],'k')
%plot([50.2 50.2],[23,43],'k')
%plot([40, 65],[37.1,37.1],'k')
%plot([40, 65],[29.9,29.9],'k')


