%% Number of static factors
function r = numstatic(x,rmax,static)

    %           _______________$$ Onatski criterion $$_____________________
if strcmp(static,'onatski')
    r                     = onatskicriterion2009(x,rmax);
    if r==0
        error('zero factors found under Onatski criterion')
    end
    
    %           _______________$$ Bai-Ng criterion $$______________________
elseif strcmp(static,'baing')
    [~, IC] = baingcriterion2002(x, rmax);
    [~,r] = min(IC(:,2));
    
    %           _______________$$ Onatski test $$__________________________
elseif strcmp(static,'onatskitest')
    r    = onatskitest(x,'s');
    
    %           _______________$$ no criterion $$__________________________
elseif isnumeric(static)
    r = static;
    
end

end
