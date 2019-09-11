function Kinetic_Correlation

% load transcription factor activity and relative metabolite abundances

warning off all

[FileName,PathName,FilterIndex]=uigetfile('*.xlsx','Select Supplementray Tables');
[status,sheets,~]= xlsfinfo([PathName FileName]);
[~,~,raw1] = xlsread([PathName FileName],6);
[~,~,raw2] = xlsread([PathName FileName],2);
TFA  = 10.^cell2mat(raw1(5:end, 2:30));  %calculate transcription factor activity from log10 transcription factor activity
Met1 = cell2mat(raw2(5:end, 5:76));

time = cell2mat(raw2(4, 5:76));     %get timepoints to average the metabolites
time2=unique(time);

for i=1:length(time2)
    id          =find(time==time2(i));
    Met(:,i)    = mean(Met1(:,id)');
end



% only timepoints where we have data for metabolites and transcription factor activities
i_met=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 20 22 24 25 26 27 28 29 30 31 32 34];
Met=Met(:,i_met);




%% go through all metabolites and transcription factor activities and calculate kinetic correlation TFA = f(Met) with and without shift

for o=1:length(Met)           
    for w=1:length(TFA)
        
        sprintf('TF: %d - Met %d', w, o)
        x=Met(o,:);
        y=TFA(w,:);
        Rsq=corr(x', y');       % calculate linear correlation to decide between activaton or inhibition kinetic

        lags=[-1 0];

        for j=1:length(lags)

            x2=[];
            y2=[];

            x2=x(1:end+lags(j));
            y2=y(-lags(j)+1:end);       % shift transcription factor activity 
            
            for ind=1:50        %run 50 fits to get best parameters for correlation

                meta_new=[];
                if Rsq>0                       % if it is positiv linear correlation use activation kinetics

                    [xData, yData] = prepareCurveData( x2, y2 );
                    ft = fittype( 'v*((x^h)/(x^h+(k^h)))', 'independent', 'x', 'dependent', 'y' );
                    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                    opts.Display = 'Off';
                    opts.MaxIter = 2000;
                    opts.Lower = [0 0 0];
                    opts.Upper = [10 Inf Inf];          % constrain hill parameter to maximum 10
                    opts.Robust = 'LAR';
                    try
                        [fitresult, gof] = fit( xData, yData, ft, opts );
                        cvs=coeffvalues(fitresult);

                        for k=1:length(x2)
                            meta_new(k)=cvs(3)*((x2(k)^cvs(1))/((x2(k)^cvs(1))+cvs(2)^cvs(1)));   % calculate the f(met) based on the kinetics
                        end
                        xi=meta_new;
                    catch
                    end
                else                            % if it is negativ linear correlation use inhibition kinetics
 
                    [xData, yData] = prepareCurveData( x2, y2 );
                    ft = fittype( 'v*((k^h)/(x^h+(k^h)))', 'independent', 'x', 'dependent', 'y' );
                    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                    opts.Display = 'Off';
                    opts.MaxIter = 2000;
                    opts.Lower =  [0 0 0];
                    opts.Upper = [10 Inf Inf];              % constrain hill parameter to maximum 10
                    opts.Robust = 'LAR';
                    try
                        [fitresult, gof] = fit( xData, yData, ft, opts );
                        cvs=coeffvalues(fitresult);

                        for k=1:length(x2)
                            meta_new(k)=cvs(3)*((cvs(2)^cvs(1))/((x2(k)^cvs(1))+cvs(2)^cvs(1)));   % calculate the f(met) based on the kinetics
                        end
                        xi=meta_new;
                    catch
                    end
                end

                if isempty(meta_new)
                else
                    Rsq_new1(ind)=corr(xi', y2');       %correlates TFA and f(Met)
                end

            end

            [mi ix]=sort(Rsq_new1,'descend');       %find best correlation of all 20 fits
            ixx=find(~isnan(mi));

            Rsq_new2(j)=Rsq_new1(ix(ixx(1)));       

        end
        
        [d ixxx]=max(Rsq_new2);     % check wether shifted or unshifted TFA gives the best result
        RSQ1(o,w)=Rsq_new2(1);      %Get correlation coefficient with TFA shiffted for one timepoint
        LAG1(o,w)=lags(1);
        RSQ0(o,w)=Rsq_new2(2);      %Get correlation coefficient without TFA shiffted for one timepoint
        LAG0(o,w)=lags(2);
        RSQ(o,w)=Rsq_new2(ixxx);    %Get correlation coefficient with TFA the best correlation
        LAG(o,w)=lags(ixxx);        %Get corrsponding shift for the best correlation
    end
    
end

save('corr_result.mat', 'RSQ1', 'RSQ0', 'RSQ', 'LAG')
