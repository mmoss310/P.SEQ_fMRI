function data=OrderYourMind_fMRI_ProcessDSVars(mtime,data)
%==========================================================================
% locProcessDSVars prepares experiment specific variables that need to be
% coded based on the preexisting values in the dataset array! Specific to
% the experiment! And, only works for dataset array! So, let's be tedious!
%==========================================================================
switch mtime
    case {'BEH'},%Behavioral file befor merging  
        % Lag 1 variables
        data.LACC=lag(data.ACC,1);
        data.LACC(data.TRIAL==1) = NaN;
        data.L1ELEMENT=lag(data.ELEMENT,1);
        
        %Lag2 variables
        data.LAG2(:,1)=zeros(size(data,1),1);
        data.LAG2(lag(data.ELEMENT,2)==data.ELEMENT,:)=1;
                
        %U1 variables
        data.UELEMENT=circshift(data.ELEMENT,[-1,1]);

        %ITEM-POS code
        data.ITEMPOSC=nan(size(data.BLOCK,1),1);
        allpairs=allcomb(unique(data.ELEMENT),unique(data.WCPOS));
        for p=1:size(allpairs,1)
            c=allpairs(p,:);
            data.ITEMPOSC(data.ELEMENT==c(1)& data.WCPOS==c(2))=p;
        end

    case {'CUTOFF'},        

        %Selection criteria for data
        %data=data(data.ACC == 1,:);
        %data=data(data.LACC == 1,:);
        data=data(data.RT > 0& data.RT < 5000,:);
        
        %Exclude subjects
        %data=data(data.SUBID~=??,:);

    case {'KGROUP'},
        %Rejoing to the original data
        kdata=v2struct(load('kscoreDataS_av.mat'));
        mergeKeys={'SUBID'};data=join(data,kdata,'key',mergeKeys);   
end


end

