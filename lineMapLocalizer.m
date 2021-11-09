classdef lineMapLocalizer < handle
	properties(Constant)
		maxErr = 0.05;
	 	minPts = 5;
	end
	properties(Access = private)

	end
	properties(Access = public)
		lines_p1 = [];
		lines_p2 = [];
		gain = 0.3;
		errThresh = 0.01;
		gradThresh = 0.0005;
    end
    
    methods
        function obj = lineMapLocalizer(lines_p1,lines_p2,gain,errThresh,gradThresh)
            obj.lines_p1 = lines_p1;
            obj.lines_p2 = lines_p2;
            obj.gain = gain;
            obj.errThresh = errThresh;
            obj.gradThresh = gradThresh;
        end
        
        function ro2 = closestSquaredDistanceToLines(obj,pi)
            pi = pi(1:2,:);
            r2Array = zeros(size(obj.lines_p1,2),size(pi,2));
            for i = 1:size(obj.lines_p1,2)
                [r2Array(i,:) , ~] = closestPointOnLineSegment(pi,...
                    obj.lines_p1(:,i),obj.lines_p2(:,i));
            end
            ro2 = min(r2Array,[],1);
        end
        
        function ids = throwOutliers(obj,pose,ptsInModelFrame)
            worldPts = pose.bToA()*ptsInModelFrame;
            r2 = obj.closestSquaredDistanceToLines(worldPts);
            ids = find(r2 > obj.maxErr*obj.maxErr);
        end
        
        function avgErr2 = fitError(obj,pose,ptsInModelFrame)
            worldPts = pose.bToA()*ptsInModelFrame;
            r2 = obj.closestSquaredDistanceToLines(worldPts);
            r2(r2 == Inf) = [];
            err2= sum(r2);
            num = length(r2);
            if(num >= lineMapLocalizer.minPts)
                avgErr2 = err2/num;
            else
                avgErr2= inf;
            end
        end
        
        function [err2_Plus0,J] = getJacobian(obj,poseVecIn,modelPts)
            err2_Plus0 = fitError(obj,poseVecIn,modelPts);
        
            eps = 1e-9;
            dpX = [eps ;0.0 ;0.0];
            newPoseVecX = poseVecIn.getPose+dpX;
            dpY = [0.0;eps ;0.0];
            newPoseVecY = poseVecIn.getPose+dpY;
            dpTh = [0.0 ;0.0 ;eps];
            newPoseVecTh = poseVecIn.getPose+dpTh;
            % Fill me in...
            delEdelX = 1/eps * (fitError(obj,newPoseVecX,modelPts) - err2_Plus0);
            delEdelY = 1/eps * (fitError(obj,newPoseVecY,modelPts) - err2_Plus0);
            delEdelTheta = 1/eps * (fitError(obj,newPoseVecTh,modelPts) - err2_Plus0);
            J = [delEdelX; delEdelY; delEdelTheta];
        end
        
        function[success, outPose] ...
                = refinePose(obj,inPose,ptsInModelFrame,maxIters)
            % refine robot pose in world (inPose) based on lidar
            % registration. Terminates if maxIters iterations is
            % exceeded or if insufficient points match the lines.
            % Even if the minimum is not found, outPose will contain 
            % anychanges that reduced the fit error. Pose changes that
            % increase fit error are not included and termination 
            % occursthereafter.
            
            % Fill me in...
            pose = inPose;
            
            success = false;
            for i=1:maxIters
                e = fitError(obj, pose, ptsInModelFrame);
                if e > obj.gradThresh
                    success = true;
                    break;
                end
                diff = -obj.gain*getJacobian(obj, pose, ptsInModelFrame);
                pose = pose + diff;
            end
            outPose = pose;
        end
        
    end
end

function [rad2, po] = closestPointOnLineSegment(pi,p1,p2)
    v1 = bsxfun(@minus,pi,p1);
    v2 = p2-p1;
    v3 = bsxfun(@minus,pi,p2);
    v1dotv2 = bsxfun(@times,v1,v2);
    v1dotv2 = sum(v1dotv2,1);
    v2dotv2 = sum(v2.*v2);
    v3dotv2 = bsxfun(@times,v3,v2);
    v3dotv2 = sum(v3dotv2,1);
    nPoints = size(pi,2);
    rad2 = zeros(1,nPoints);
    po = zeros(2,nPoints);
    % Closest is on segment
    flag1 = v1dotv2 > 0.0 & v3dotv2 < 0.0;
    if any(flag1)
        scale = v1dotv2/v2dotv2;
        temp = bsxfun(@plus,v2*scale,[p1(1) ; p1(2)]);
        po(:,flag1) = temp(:,flag1);
        dx = pi(1,flag1)-po(1,flag1);
        dy = pi(2,flag1)-po(2,flag1);
        rad2(flag1) = dx.*dx+dy.*dy;
    end% Closest is first endpoint
    flag2 = v1dotv2 <= 0.0;
    if any(flag2)
        temp = bsxfun(@times,ones(2,sum(flag2)),[p1(1); p1(2)]);
        po(:,flag2) = temp;
        rad2(flag2) = inf;
    end% Closest is second endpoint
    flag3 = ~flag1 & ~flag2;
    if any(flag3)
        temp = bsxfun(@times,ones(2,sum(flag3)),[p2(1); p2(2)]);
        po(:,flag3) = temp;
        rad2(flag3) = inf;
    end
end