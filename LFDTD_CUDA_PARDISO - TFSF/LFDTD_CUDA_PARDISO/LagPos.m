classdef LagPos
    
    methods (Static)
        
        function [im] = Plot_E_Field(plotFieldEx,plotFieldEy,plotFieldEz,plotSec,recordE,nodeNum,nx,ny,nz,plotStep,plotLayer)
            
            exPlot = zeros(nx,ny+1,nz+1,plotStep);
            eyPlot = zeros(nx+1,ny,nz+1,plotStep);
            ezPlot = zeros(nx+1,ny+1,nz,plotStep);
            
            % Ex Field
            
            for i=1:nx
                for j=1:ny+1
                    for k=1:nz+1
                        
                        ex111 = nodeNum(i,j,k,1);
                        exPlot(i,j,k,:) = recordE(ex111,:);
                    end
                end
            end
            
            % Ey Field
            
            for i=1:nx+1
                for j=1:ny
                    for k=1:nz+1
                        
                        ey111 = nodeNum(i,j,k,2);
                        eyPlot(i,j,k,:) = recordE(ey111,:);
                    end
                end
            end
            
            % Ez Field
            
            for i=1:nx+1
                for j=1:ny+1
                    for k=1:nz
                        
                        ez111 = nodeNum(i,j,k,3);
                        ezPlot(i,j,k,:) = recordE(ez111,:);
                    end
                end
            end
            
            % Plot Ex
            
            if plotFieldEx==1
                
                if plotSec==1
                    
                    cmax = max(reshape(exPlot(plotLayer,:,:,:),(ny+1)*(nz+1)*plotStep,1)); % Max E field value
                    cmin = min(reshape(exPlot(plotLayer,:,:,:),(ny+1)*(nz+1)*plotStep,1)); % Min E field value
                    plotLimit = max([abs(cmin),abs(cmax)]);
                    figure;
                    F(plotStep) = struct('cdata',[],'colormap',[]);
                    for i=1:plotStep
                        surf((reshape(exPlot(plotLayer,:,:,i),ny+1,nz+1))');
                        shading flat;
                        clim([-plotLimit,plotLimit]); % Value range of the colorbar
                        colorbar;
                        axis([1,ny+1,1,nz+1,-plotLimit,plotLimit]);
                        title(strcat('Ex; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                        F(i) = getframe;
                    end
                    
                elseif plotSec==2
                    
                    cmax = max(reshape(exPlot(:,plotLayer,:,:),nx*(nz+1)*plotStep,1)); % Max E field value
                    cmin = min(reshape(exPlot(:,plotLayer,:,:),nx*(nz+1)*plotStep,1)); % Min E field value
                    plotLimit = max([abs(cmin),abs(cmax)]);
                    figure;
                    F(plotStep) = struct('cdata',[],'colormap',[]);
                    for i=1:plotStep
                        surf((reshape(exPlot(:,plotLayer,:,i),nx,nz+1))');
                        shading flat;
                        clim([-plotLimit,plotLimit]); % Value range of the colorbar
                        colorbar;
                        axis([1,nx,1,nz+1,-plotLimit,plotLimit]);
                        title(strcat('Ex; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                        F(i) = getframe;
                    end
                    
                elseif plotSec==3
                    
                    cmax = max(reshape(exPlot(:,:,plotLayer,:),nx*(ny+1)*plotStep,1)); % Max E field value
                    cmin = min(reshape(exPlot(:,:,plotLayer,:),nx*(ny+1)*plotStep,1)); % Min E field value
                    plotLimit = max([abs(cmin),abs(cmax)]);
                    figure;
                    F(plotStep) = struct('cdata',[],'colormap',[]);
                    im = cell(plotStep,1);
                    for i=1:plotStep
                        surf((reshape(exPlot(:,:,plotLayer,i),nx,ny+1))');
                        shading flat;
                        clim([-plotLimit,plotLimit]); % Value range of the colorbar
                        colorbar;
                        axis([1,nx,1,ny+1,-plotLimit,plotLimit]);
                        title(strcat('Ex; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                        view(3);
                        drawnow;
                    end
                    
                else
                    
                    fprintf('Plot section setting error!\n');
                    return;
                    
                end
                
            end
            
            % Plot Ey
            
            if plotFieldEy==1
                
                if plotSec==1
                    
                    cmax = max(reshape(eyPlot(plotLayer,:,:,:),ny*(nz+1)*plotStep,1)); % Max E field value
                    cmin = min(reshape(eyPlot(plotLayer,:,:,:),ny*(nz+1)*plotStep,1)); % Min E field value
                    plotLimit = max([abs(cmin),abs(cmax)]);
                    figure;
                    F(plotStep) = struct('cdata',[],'colormap',[]);
                    for i=1:plotStep
                        surf((reshape(eyPlot(plotLayer,:,:,i),ny,nz+1))');
                        shading flat;
                        clim([-plotLimit,plotLimit]); % Value range of the colorbar
                        colorbar;
                        axis([1,ny,1,nz+1,-plotLimit,plotLimit]);
                        title(strcat('Ey; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                        F(i) = getframe;
                    end
                    
                elseif plotSec==2
                    
                    cmax = max(reshape(eyPlot(:,plotLayer,:,:),(nx+1)*(nz+1)*plotStep,1)); % Max E field value
                    cmin = min(reshape(eyPlot(:,plotLayer,:,:),(nx+1)*(nz+1)*plotStep,1)); % Min E field value
                    plotLimit = max([abs(cmin),abs(cmax)]);
                    figure;
                    F(plotStep) = struct('cdata',[],'colormap',[]);
                    for i=1:plotStep
                        surf((reshape(eyPlot(:,plotLayer,:,i),nx+1,nz+1))');
                        shading flat;
                        clim([-plotLimit,plotLimit]); % Value range of the colorbar
                        colorbar;
                        axis([1,nx+1,1,nz+1,-plotLimit,plotLimit]);
                        title(strcat('Ey; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                        F(i) = getframe;
                    end
                    
                elseif plotSec==3
                    
                    cmax = max(reshape(eyPlot(:,:,plotLayer,:),(nx+1)*ny*plotStep,1)); % Max E field value
                    cmin = min(reshape(eyPlot(:,:,plotLayer,:),(nx+1)*ny*plotStep,1)); % Min E field value
                    plotLimit = max([abs(cmin),abs(cmax)]);
                    fig = figure;
                    F(plotStep) = struct('cdata',[],'colormap',[]);
                    im = cell(plotStep,1);
                    for i=1:plotStep
                        surf((reshape(eyPlot(:,:,plotLayer,i),nx+1,ny))');
                        shading flat;
                        clim([-plotLimit,plotLimit]); % Value range of the colorbar
                        colorbar;
                        axis([1,nx+1,1,ny,-plotLimit,plotLimit]);
                        title(strcat('Ey; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                        drawnow;
                        frame = getframe(fig);
                        im{i} = frame2im(frame);
                        F(i) = getframe(gca);
                    end
                    
                else
                    
                    fprintf('Plot section setting error!\n');
                    return;
                    
                end
                
            end
            
            % Plot Ez
            
            if plotFieldEz==1
                
                if plotSec==1
                    
                    cmax = max(reshape(ezPlot(plotLayer,:,:,:),(ny+1)*nz*plotStep,1)); % Max E field value
                    cmin = min(reshape(ezPlot(plotLayer,:,:,:),(ny+1)*nz*plotStep,1)); % Min E field value
                    plotLimit = max([abs(cmin),abs(cmax)]);
                    fig = figure;
                    F(plotStep) = struct('cdata',[],'colormap',[]);
                    im = cell(plotStep,1);
                    for i=1:plotStep
                        surf((reshape(ezPlot(plotLayer,:,:,i),ny+1,nz))');
                        shading flat;
                        clim([-plotLimit,plotLimit]); % Value range of the colorbar
                        colorbar;
                        axis([1,ny+1,1,nz,-plotLimit,plotLimit]);
                        title(strcat('Ez; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                        frame = getframe(fig);
                        im{i} = frame2im(frame);
                        F(i) = getframe(gca);
                    end
                    
                elseif plotSec==2
                    
                    cmax = max(reshape(ezPlot(:,plotLayer,:,:),(nx+1)*nz*plotStep,1)); % Max E field value
                    cmin = min(reshape(ezPlot(:,plotLayer,:,:),(nx+1)*nz*plotStep,1)); % Min E field value
                    plotLimit = max([abs(cmin),abs(cmax)]);
                    figure;
                    F(plotStep) = struct('cdata',[],'colormap',[]);
                    for i=1:plotStep
                        surf((reshape(ezPlot(:,plotLayer,:,i),nx+1,nz))');
                        shading flat;
                        clim([-plotLimit,plotLimit]); % Value range of the colorbar
                        colorbar;
                        axis([1,nx+1,1,nz,-plotLimit,plotLimit]);
                        title(strcat('Ez; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                        F(i) = getframe;
                    end
                    
                elseif plotSec==3
                    
                    cmax = max(reshape(ezPlot(:,:,plotLayer,:),(nx+1)*(ny+1)*plotStep,1)); % Max E field value
                    cmin = min(reshape(ezPlot(:,:,plotLayer,:),(nx+1)*(ny+1)*plotStep,1)); % Min E field value
                    plotLimit = max([abs(cmin),abs(cmax)]);
                    fig = figure;
                    F(plotStep) = struct('cdata',[],'colormap',[]);
                    im = cell(plotStep,1);
                    for i=1:plotStep
                        surf((reshape(ezPlot(:,:,plotLayer,i),nx+1,ny+1))');
                        shading flat;
                        clim([-plotLimit,plotLimit]); % Value range of the colorbar
                        colorbar;
                        axis([1,nx+1,1,ny+1,-plotLimit,plotLimit]);
                        title(strcat('Ez; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                        frame = getframe(fig);
                        im{i} = frame2im(frame);
                        F(i) = getframe(gca);
                    end
                    
                else
                    
                    fprintf('Plot section setting error!\n');
                    return;
                    
                end
        
            end
        end
        
        function [im] = Plot_H_Field(plotFieldHx,plotFieldHy,plotFieldHz,plotSec,recordHx,recordHy,recordHz,nx,ny,nz,plotStep,plotLayer)
            
                % Plot Hx
    
                if plotFieldHx == 1
                    
                    if plotSec == 1
                        cmax = max(reshape(recordHx(plotLayer,:,:,:),(ny)*(nz)*plotStep,1)); % Max H field value
                        cmin = min(reshape(recordHx(plotLayer,:,:,:),(ny)*(nz)*plotStep,1)); % Min H field value
                        plotLimit = max([abs(cmin),abs(cmax)]);
                        figure;
                        F(plotStep) = struct('cdata',[],'colormap',[]);
                        for i=1:plotStep
                            surf((reshape(recordHx(plotLayer,:,:,i),ny,nz))');
                            shading flat;
                            clim([-plotLimit,plotLimit]); % Value range of the colorbar
                            colorbar;
                            axis([1,ny+1,1,nz+1,-plotLimit,plotLimit]);
                            title(strcat('Hx; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                            F(i) = getframe;
                        end
                        
                    elseif plotSec==2
                        
                        cmax = max(reshape(recordHx(:,plotLayer,:,:),nx*(nz)*plotStep,1)); % Max E field value
                        cmin = min(reshape(recordHx(:,plotLayer,:,:),nx*(nz)*plotStep,1)); % Min E field value
                        plotLimit = max([abs(cmin),abs(cmax)]);
                        figure;
                        F(plotStep) = struct('cdata',[],'colormap',[]);
                        for i=1:plotStep
                            surf((reshape(recordHx(:,plotLayer,:,i),nx,nz))');
                            shading flat;
                            clim([-plotLimit,plotLimit]); % Value range of the colorbar
                            colorbar;
                            axis([1,nx,1,nz+1,-plotLimit,plotLimit]);
                            title(strcat('Hx; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                            F(i) = getframe;
                        end
                        
                    elseif plotSec==3
                        
                        cmax = max(reshape(recordHx(:,:,plotLayer,:),nx*(ny)*plotStep,1)); % Max E field value
                        cmin = min(reshape(recordHx(:,:,plotLayer,:),nx*(ny)*plotStep,1)); % Min E field value
                        plotLimit = max([abs(cmin),abs(cmax)]);
                        figure;
                        F(plotStep) = struct('cdata',[],'colormap',[]);
                        for i=1:plotStep
                            surf((reshape(recordHx(:,:,plotLayer,i),nx,ny))');
                            shading flat;
                            clim([-plotLimit,plotLimit]); % Value range of the colorbar
                            colorbar;
                            axis([1,nx+1,1,ny+1,-plotLimit,plotLimit]);
                            title(strcat('Hx; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                            F(i) = getframe;
                        end
                    end
                    
                end
                
                % Plot Hy
                if plotFieldHy == 1
                    
                    if plotSec == 1
                        cmax = max(reshape(recordHy(plotLayer,:,:,:),ny*(nz)*plotStep,1)); % Max E field value
                        cmin = min(reshape(recordHy(plotLayer,:,:,:),ny*(nz)*plotStep,1)); % Min E field value
                        plotLimit = max([abs(cmin),abs(cmax)]);
                        figure;
                        F(plotStep) = struct('cdata',[],'colormap',[]);
                        for i=1:plotStep
                            surf((reshape(recordHy(plotLayer,:,:,i),ny,nz))');
                            shading flat;
                            clim([-plotLimit,plotLimit]); % Value range of the colorbar
                            colorbar;
                            axis([1,ny,1,nz+1,-plotLimit,plotLimit]);
                            title(strcat('Hy; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                            F(i) = getframe;
                        end
                        
                    elseif plotSec==2
                        
                        cmax = max(reshape(recordHy(:,plotLayer,:,:),(nx)*(nz)*plotStep,1)); % Max E field value
                        cmin = min(reshape(recordHy(:,plotLayer,:,:),(nx)*(nz)*plotStep,1)); % Min E field value
                        plotLimit = max([abs(cmin),abs(cmax)]);
                        figure;
                        F(plotStep) = struct('cdata',[],'colormap',[]);
                        for i=1:plotStep
                            surf((reshape(recordHy(:,plotLayer,:,i),nx,nz))');
                            shading flat;
                            clim([-plotLimit,plotLimit]); % Value range of the colorbar
                            colorbar;
                            axis([1,nx+1,1,nz+1,-plotLimit,plotLimit]);
                            title(strcat('Hy; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                            F(i) = getframe;
                        end
                        
                        
                    elseif plotSec==3
                        
                        cmax = max(reshape(recordHy(:,:,plotLayer,:),(nx)*ny*plotStep,1)); % Max E field value
                        cmin = min(reshape(recordHy(:,:,plotLayer,:),(nx)*ny*plotStep,1)); % Min E field value
                        plotLimit = max([abs(cmin),abs(cmax)]);
                        fig = figure;
                        F(plotStep) = struct('cdata',[],'colormap',[]);
                        im = cell(plotStep,1);
                        for i=1:plotStep
                            surf((reshape(recordHy(:,:,plotLayer,i),nx,ny))');
                            shading flat;
                            clim([-plotLimit,plotLimit]); % Value range of the colorbar
                            colorbar;
                            axis([1,nx+1,1,ny,-plotLimit,plotLimit]);
                            title(strcat('Hy; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                            F(i) = getframe;
                            drawnow;
                            frame = getframe(fig);
                            im{i} = frame2im(frame);
                            F(i) = getframe(gca);
                        end
                    end
                    
                    
                end
                
                % Plot Hz
                
                if plotFieldHz == 1
                    
                    if plotSec == 1
                        
                        cmax = max(reshape(recordHz(plotLayer,:,:,:),(ny)*nz*plotStep,1)); % Max E field value
                        cmin = min(reshape(recordHz(plotLayer,:,:,:),(ny)*nz*plotStep,1)); % Min E field value
                        plotLimit = max([abs(cmin),abs(cmax)]);
                        figure;
                        F(plotStep) = struct('cdata',[],'colormap',[]);
                        for i=1:plotStep
                            surf((reshape(recordHz(plotLayer,:,:,i),ny,nz))');
                            shading flat;
                            clim([-plotLimit,plotLimit]); % Value range of the colorbar
                            colorbar;
                            axis([1,ny+1,1,nz,-plotLimit,plotLimit]);
                            title(strcat('Hz; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                            F(i) = getframe;
                        end
                        
                    elseif plotSec==2
                        
                        cmax = max(reshape(recordHz(:,plotLayer,:,:),(nx)*nz*plotStep,1)); % Max E field value
                        cmin = min(reshape(recordHz(:,plotLayer,:,:),(nx)*nz*plotStep,1)); % Min E field value
                        plotLimit = max([abs(cmin),abs(cmax)]);
                        figure;
                        F(plotStep) = struct('cdata',[],'colormap',[]);
                        for i=1:plotStep
                            surf((reshape(recordHz(:,plotLayer,:,i),nx,nz))');
                            shading flat;
                            clim([-plotLimit,plotLimit]); % Value range of the colorbar
                            colorbar;
                            axis([1,nx+1,1,nz,-plotLimit,plotLimit]);
                            title(strcat('Hz; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                            F(i) = getframe;
                        end
                        
                    elseif plotSec==3
                        
                        cmax = max(reshape(recordHz(:,:,plotLayer,:),(nx)*(ny)*plotStep,1)); % Max E field value
                        cmin = min(reshape(recordHz(:,:,plotLayer,:),(nx)*(ny)*plotStep,1)); % Min E field value
                        plotLimit = max([abs(cmin),abs(cmax)]);
                        figure;
                        F(plotStep) = struct('cdata',[],'colormap',[]);
                        for i=1:plotStep
                            surf((reshape(recordHz(:,:,plotLayer,i),nx,ny))');
                            shading flat;
                            clim([-plotLimit,plotLimit]); % Value range of the colorbar
                            colorbar;
                            axis([1,nx+1,1,ny+1,-plotLimit,plotLimit]);
                            title(strcat('Hz; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                            F(i) = getframe;
                            %                 view(0,0);
                            %                 drawnow;
                        end
                        
                    end
                end
        end
        
        function Plot_1D_Field(recordHy_1D,nnode_1D,plotStep)
            
            cmax = max(reshape(recordHy_1D,numel(recordHy_1D),1)); % Max E field value
            cmin = min(reshape(recordHy_1D,numel(recordHy_1D),1)); % Min E field Value
            plotLimit = max([abs(cmin),abs(cmax)]);
            figure;
            for i=1:plotStep
                plot(reshape(recordHy_1D(1,1,:,i),[],1),'LineWidth',2);
                axis([1,nnode_1D-1,-plotLimit,plotLimit]);
                grid on;
                grid minor;
                title(strcat('Hy 1D; ',';Total Plot Step = ',num2str(plotStep),';Plot Step = ',num2str(i)));
                drawnow;
            end
        end
        
        function Plot_E_Waveform(tArray,fac,dt,tStep,plotStep,recordE,nodeNum,nx,ny,nz,nx_TFSF,ny_TFSF,nz_TFSF)
            figure;
            
            % Time Domain Plot
            
            fieldplot = readtable("../CST Files/Freq_CST_Cylinder_Modulated_scatter.txt");
            xpos = fieldplot{:,1};
            ypos = fieldplot{:,2};
            fixpos = (xpos == -171) & (ypos == 0);
            Exfield = fieldplot{:,4};
            Ex_plot = Exfield(fixpos);
            timeArray = linspace(0,9.9,length(Ex_plot));
            plot(timeArray,Ex_plot,'Linewidth',2);
            hold on;
%             recordE(nodeNum(nx-nx_TFSF-1,ny/2,floor(nz/2),1),342:end) = 0.5*recordE(nodeNum(nx-nx_TFSF-1,ny/2,floor(nz/2),1),342:end);
            % recordE(nodeNum(nx-nx_TFSF,ny/2,floor(nz/2),1),386:end) = 0.5*recordE(nodeNum(nx-nx_TFSF,ny/2,floor(nz/2),1),386:end);
            plot(tArray(1:fac:tStep)/10^(-9),recordE(nodeNum(nx-nx_TFSF,ny/2,floor(nz/2),1),:),'Linewidth',2);
            legend('CST','LFDTD');
            grid on;
            grid minor;
            xlabel('Time (ns)');
            ylabel('E (V/m)');
            axis tight;
            title("E_x Field Plot Time Domain");
        end
        
        function [im] = Plot_VecE_Field(plotSec,recordE,nodeNum,nx,ny,nz,plotStep,plotLayer)
            
            exPlot = zeros(nx,ny+1,nz+1,plotStep);
            eyPlot = zeros(nx+1,ny,nz+1,plotStep);
            ezPlot = zeros(nx+1,ny+1,nz,plotStep);
            
            EX = zeros(ny+1,nx,plotStep);
            EY = zeros(ny,nx+1,plotStep);
            EZ = zeros(ny+1,nx+1,plotStep);
            
            % Ex Field
            
            for i=1:nx
                for j=1:ny+1
                    for k=1:nz+1
                        
                        ex111 = nodeNum(i,j,k,1);
                        exPlot(i,j,k,:) = recordE(ex111,:);
                    end
                end
            end
            
            % Ey Field
            
            for i=1:nx+1
                for j=1:ny
                    for k=1:nz+1
                        
                        ey111 = nodeNum(i,j,k,2);
                        eyPlot(i,j,k,:) = recordE(ey111,:);
                    end
                end
            end
            
            % Ez Field
            
            for i=1:nx+1
                for j=1:ny+1
                    for k=1:nz
                        
                        ez111 = nodeNum(i,j,k,3);
                        ezPlot(i,j,k,:) = recordE(ez111,:);
                    end
                end
            end

            if plotSec == 1
                cmax = max(reshape(exPlot(plotLayer,:,:,:),(ny+1)*(nz+1)*plotStep,1)); % Max E field value
                cmin = min(reshape(exPlot(plotLayer,:,:,:),(ny+1)*(nz+1)*plotStep,1)); % Min E field value
                plotLimit = max([abs(cmin),abs(cmax)]);
                fig = figure;
                F(plotStep) = struct('cdata',[],'colormap',[]);
                im = cell(plotStep,1);
                for i = 1:plotStep
                    EX_x = (reshape(exPlot(plotLayer,:,:,i),ny+1,nz+1))';
                    EY_y = (reshape(eyPlot(plotLayer,:,:,i),ny,nz+1))';
                    EZ_z = (reshape(ezPlot(plotLayer,:,:,i),ny+1,nz))';
                    quiver3(zeros(nz,ny),EY_y(1:nz,1:ny),EZ_z(:,1:ny),EX_x(1:nz,1:ny),-2); % set Ex in upward direction
                    axis([0,nx+1,0,nz,-plotLimit,plotLimit]);
                    title(strcat('Ez; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                    view(3);
                    drawnow;
                    frame = getframe(fig);
                    im{i} = frame2im(frame);
                    F(i) = getframe(gca);
                end
            end
            
            if plotSec == 2
                cmax = max(reshape(ezPlot(:,plotLayer,:,:),(nx+1)*(nz)*plotStep,1)); % Max E field value
                cmin = min(reshape(ezPlot(:,plotLayer,:,:),(nx+1)*(nz)*plotStep,1)); % Min E field value
                plotLimit = max([abs(cmin),abs(cmax)]);
                fig = figure;
                F(plotStep) = struct('cdata',[],'colormap',[]);
                im = cell(plotStep,1);
                for i = 1:plotStep
                    EX_x = (reshape(exPlot(:,plotLayer,:,i),nx,nz+1))';
                    EY_y = (reshape(eyPlot(:,plotLayer,:,i),nx+1,nz+1))';
                    EZ_z = (reshape(ezPlot(:,plotLayer,:,i),nx+1,nz))';
                    quiver3(zeros(nz,nx),EX_x(1:nz,1:nx),EY_y(1:nz,1:nx),EZ_z(:,1:nx),-2); % set Ey in upward direction
                    axis([0,nx+1,0,nz,-plotLimit,plotLimit]);
                    title(strcat('Ez; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                    view(3);
                    drawnow;
                    frame = getframe(fig);
                    im{i} = frame2im(frame);
                    F(i) = getframe(gca);
                end
            end
            
            if plotSec == 3
                cmax = max(reshape(exPlot(:,:,plotLayer,:),(nx)*(ny+1)*plotStep,1)); % Max E field value
                cmin = min(reshape(exPlot(:,:,plotLayer,:),(nx)*(ny+1)*plotStep,1)); % Min E field value
                plotLimit = max([abs(cmin),abs(cmax)]);
                fig = figure;
                F(plotStep) = struct('cdata',[],'colormap',[]);
                im = cell(plotStep,1);
                for i = 1:plotStep
                    EX(:,:,i) = (reshape(exPlot(:,:,plotLayer,i),nx,ny+1))';
                    EY(:,:,i) = (reshape(eyPlot(:,:,plotLayer,i),nx+1,ny))';
                    EZ(:,:,i) = (reshape(ezPlot(:,:,plotLayer,i),nx+1,ny+1))';
                    quiver3(zeros(ny,nx),EX(1:ny,:,i),EY(:,1:nx,i),EZ(1:ny,1:nx,i),-1);
                    axis([0,nx+1,0,ny+1,-plotLimit,plotLimit]);
                    title(strcat('Ez; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                    view(3);
                    drawnow;
                    frame = getframe(fig);
                    im{i} = frame2im(frame);
                    F(i) = getframe(gca);
                end
            end
        end
        
        function [im] = Plot_VecH_Field(plotSec,plotLayer,nx,ny,nz,plotStep,recordHx,recordHy,recordHz)
            
            HX = zeros(nx,ny,plotStep);
            HY = zeros(nx,ny,plotStep);
            HZ = zeros(nx,ny,plotStep);
            
            if plotSec == 3
                
                cmax = max(reshape(recordHy(:,:,plotLayer,:),(nx)*(ny)*plotStep,1)); % Max H field value
                cmin = min(reshape(recordHy(:,:,plotLayer,:),(nx)*(ny)*plotStep,1)); % Min H field value
                plotLimit = max([abs(cmin),abs(cmax)]);
                fig = figure;
                F(plotStep) = struct('cdata',[],'colormap',[]);
                im = cell(plotStep,1);
                for i = 1:plotStep
                    HX(:,:,i) = (reshape(recordHx(:,:,plotLayer,i),nx,ny))';
                    HY(:,:,i) = (reshape(recordHy(:,:,plotLayer,i),nx,ny))';
                    HZ(:,:,i) = (reshape(recordHz(:,:,plotLayer,i),nx,ny))';
                    quiver3(zeros(ny,nx),HX(:,:,i),HY(:,:,i),HZ(:,:,i),1);
                    axis([0,nx+1,0,ny+1,-plotLimit,plotLimit]);
                    title(strcat('Hz; ','Plot Section = ',num2str(plotSec),'; Plot Layer = ',num2str(plotLayer),'; Total Plot Step = ',num2str(plotStep),'; Plot Step = ',num2str(i)));
                    view(3);
                    drawnow;
                    frame = getframe(fig);
                    im{i} = frame2im(frame);
                    F(i) = getframe(gca);
                end
            end
        end
        
        function Plot_Mesh(nx,ny,nz,lx,ly,lz,dx,dy,dz,eps,sigma,eps0)
            
            figure;
            set(gcf,'Position',[50,150,1500,800], 'color','w');
            meshStyle = 3; % 1-All vacuum cells without edge; 2-All vacuum with edge; 3-Boundary vacuum cell with edge
%             vacColor = 'w'; % Color of vacuum cells
            
            % Mesh
            hold on;
            
            for i=1:nx
                for j=1:ny
                    for k=8:8
                        
                        if i==1
                            x1 = 0;
                        else
                            x1 = sum(dx(1:i-1));
                        end
                        x2 = x1+dx(i);
                        
                        x4 = x1;
                        x5 = x1;
                        x8 = x1;
                        
                        x3 = x2;
                        x6 = x2;
                        x7 = x2;
                        
                        
                        if j==1
                            y1 = 0;
                        else
                            y1 = sum(dy(1:j-1));
                        end
                        y4 = y1+dy(j);
                        
                        y2 = y1;
                        y5 = y1;
                        y6 = y1;
                        
                        y3 = y4;
                        y7 = y4;
                        y8 = y4;
                        
                        if k==1
                            z1 = 0;
                        else
                            z1 = sum(dz(1:k-1));
                        end
                        z5 = z1+dz(k);
                        
                        z2 = z1;
                        z3 = z1;
                        z4 = z1;
                        
                        z6 = z5;
                        z7 = z5;
                        z8 = z5;
                        
                        if eps(i,j,k)~=eps0 || sigma(i,j,k)~=0
                            
                            h1 = patch([x1,x2,x3,x4],[y1,y2,y3,y4],[z1,z2,z3,z4],'y','LineWidth',0.1);
                            h2 = patch([x1,x5,x6,x2],[y1,y5,y6,y2],[z1,z5,z6,z2],'y','LineWidth',0.1);
                            h3 = patch([x6,x2,x3,x7],[y6,y2,y3,y7],[z6,z2,z3,z7],'y','LineWidth',0.1);
                            h4 = patch([x7,x3,x4,x8],[y7,y3,y4,y8],[z7,z3,z4,z8],'y','LineWidth',0.1);
                            h5 = patch([x5,x1,x4,x8],[y5,y1,y4,y8],[z5,z1,z4,z8],'y','LineWidth',0.1);
                            h6 = patch([x5,x6,x7,x8],[y5,y6,y7,y8],[z5,z6,z7,z8],'y','LineWidth',0.1);
                            
                            alpha(h1,0.9);
                            alpha(h2,0.9);
                            alpha(h3,0.9);
                            alpha(h4,0.9);
                            alpha(h5,0.9);
                            alpha(h6,0.9);
                            
                        else
                            
                            if meshStyle==1
                                
                                h1 = patch([x1,x2,x3,x4],[y1,y2,y3,y4],[z1,z2,z3,z4],'r','LineWidth',0.1,'EdgeColor','none');
                                h2 = patch([x1,x5,x6,x2],[y1,y5,y6,y2],[z1,z5,z6,z2],'r','LineWidth',0.1,'EdgeColor','none');
                                h3 = patch([x6,x2,x3,x7],[y6,y2,y3,y7],[z6,z2,z3,z7],'r','LineWidth',0.1,'EdgeColor','none');
                                h4 = patch([x7,x3,x4,x8],[y7,y3,y4,y8],[z7,z3,z4,z8],'r','LineWidth',0.1,'EdgeColor','none');
                                h5 = patch([x5,x1,x4,x8],[y5,y1,y4,y8],[z5,z1,z4,z8],'r','LineWidth',0.1,'EdgeColor','none');
                                h6 = patch([x5,x6,x7,x8],[y5,y6,y7,y8],[z5,z6,z7,z8],'r','LineWidth',0.1,'EdgeColor','none');
                                
                            elseif meshStyle==2
                                
                                h1 = patch([x1,x2,x3,x4],[y1,y2,y3,y4],[z1,z2,z3,z4],'r','LineWidth',0.1,'EdgeColor','k');
                                h2 = patch([x1,x5,x6,x2],[y1,y5,y6,y2],[z1,z5,z6,z2],'r','LineWidth',0.1,'EdgeColor','k');
                                h3 = patch([x6,x2,x3,x7],[y6,y2,y3,y7],[z6,z2,z3,z7],'r','LineWidth',0.1,'EdgeColor','k');     
                                h4 = patch([x7,x3,x4,x8],[y7,y3,y4,y8],[z7,z3,z4,z8],'r','LineWidth',0.1,'EdgeColor','k');
                                h5 = patch([x5,x1,x4,x8],[y5,y1,y4,y8],[z5,z1,z4,z8],'r','LineWidth',0.1,'EdgeColor','k');
                                h6 = patch([x5,x6,x7,x8],[y5,y6,y7,y8],[z5,z6,z7,z8],'r','LineWidth',0.1,'EdgeColor','k');
                                
                            elseif meshStyle==3
                                
                                if x1==0
                                    h5 = patch([x5,x1,x4,x8],[y5,y1,y4,y8],[z5,z1,z4,z8],'w','LineWidth',0.1,'EdgeColor','k');
                                    h3 = patch([x6,x2,x3,x7],[y6,y2,y3,y7],[z6,z2,z3,z7],'w','LineWidth',0.1,'EdgeColor','none');
                                elseif x2==lx
                                    h5 = patch([x5,x1,x4,x8],[y5,y1,y4,y8],[z5,z1,z4,z8],'w','LineWidth',0.1,'EdgeColor','none');
                                    h3 = patch([x6,x2,x3,x7],[y6,y2,y3,y7],[z6,z2,z3,z7],'w','LineWidth',0.1,'EdgeColor','k');
                                else
                                    h5 = patch([x5,x1,x4,x8],[y5,y1,y4,y8],[z5,z1,z4,z8],'w','LineWidth',0.1,'EdgeColor','none');
                                    h3 = patch([x6,x2,x3,x7],[y6,y2,y3,y7],[z6,z2,z3,z7],'w','LineWidth',0.1,'EdgeColor','none');
                                end
                                
                                if y1==0
                                    h2 = patch([x1,x5,x6,x2],[y1,y5,y6,y2],[z1,z5,z6,z2],'w','LineWidth',0.1,'EdgeColor','k');
                                    h4 = patch([x7,x3,x4,x8],[y7,y3,y4,y8],[z7,z3,z4,z8],'w','LineWidth',0.1,'EdgeColor','none');
                                elseif y4==ly
                                    h2 = patch([x1,x5,x6,x2],[y1,y5,y6,y2],[z1,z5,z6,z2],'w','LineWidth',0.1,'EdgeColor','none');
                                    h4 = patch([x7,x3,x4,x8],[y7,y3,y4,y8],[z7,z3,z4,z8],'w','LineWidth',0.1,'EdgeColor','k');
                                else
                                    h2 = patch([x1,x5,x6,x2],[y1,y5,y6,y2],[z1,z5,z6,z2],'w','LineWidth',0.1,'EdgeColor','none');
                                    h4 = patch([x7,x3,x4,x8],[y7,y3,y4,y8],[z7,z3,z4,z8],'w','LineWidth',0.1,'EdgeColor','none');
                                end
                                
                                if z1==0
                                    h1 = patch([x1,x2,x3,x4],[y1,y2,y3,y4],[z1,z2,z3,z4],'w','LineWidth',0.1,'EdgeColor','k');
                                    h6 = patch([x5,x6,x7,x8],[y5,y6,y7,y8],[z5,z6,z7,z8],'w','LineWidth',0.1,'EdgeColor','none');
                                elseif z5==lz
                                    h1 = patch([x1,x2,x3,x4],[y1,y2,y3,y4],[z1,z2,z3,z4],'w','LineWidth',0.1,'EdgeColor','none');
                                    h6 = patch([x5,x6,x7,x8],[y5,y6,y7,y8],[z5,z6,z7,z8],'w','LineWidth',0.1,'EdgeColor','k');
                                else
                                    h1 = patch([x1,x2,x3,x4],[y1,y2,y3,y4],[z1,z2,z3,z4],'w','LineWidth',0.1,'EdgeColor','none');
                                    h6 = patch([x5,x6,x7,x8],[y5,y6,y7,y8],[z5,z6,z7,z8],'w','LineWidth',0.1,'EdgeColor','none');
                                end
                                
                            else
                                
                                fprintf('Mesh style setting error!\n');
                                return;
                                
                            end
                            
                            alpha(h1,0.03);
                            alpha(h2,0.03);
                            alpha(h3,0.03);
                            alpha(h4,0.03);
                            alpha(h5,0.03);
                            alpha(h6,0.03);
                            
                        end
                        
                    end
                end
            end
            
            xlabel('x (m)');
            ylabel('y (m)');
            zlabel('z (m)');
            title('Mesh ');
            axis equal;
            axis([0,sum(dx(1:nx)),0,sum(dy(1:ny)),0,sum(dz(1:nz))]);
            %     view(-37.5,30);
            view(0,90);
        end
        
    end
end