function [fx,gx] = AutoDiff(f,x,reldelta,method)
fx = f(x);
gx = zeros(length(fx),length(x));
switch method
    case 'stencil'
        for i=1:numel(x)
            h = max(reldelta*abs(x(i)),1e-6);
            
            newx = x;
            newx(i) = newx(i) + h;
            ffh = f(newx);
            
            newx = x;
            newx(i) = newx(i) + 2*h;
            ff2h = f(newx);
            
            newx = x;
            newx(i) = newx(i) - h;
            fbh = f(newx);
            
            newx = x;
            newx(i) = newx(i) - 2*h;
            fb2h = f(newx);
            
            gx(:,i) = (-ff2h + 8*ffh - 8*fbh + fb2h) / (12*h);
        end
    case 'central'
        for i=1:numel(x)
            delta = max(reldelta*abs(x(i)),1e-6)*0.5;
            
            % forward difference
            forwardx = x;
            forwardx(i) = x(i)+delta;
            forwardfx = f(forwardx);
            % backward difference
            backwardx = x;
            backwardx(i) = x(i)-delta;
            backwardfx = f(backwardx);
            gx(:,i) = (forwardfx-backwardfx)/(2*delta);
        end
    case 'forward'
        for i=1:numel(x)
            delta = max(reldelta*abs(x(i)),1e-6);
            
            % forward difference
            forwardx = x;
            forwardx(i) = x(i)+delta;
            forwardfx = f(forwardx);
            gx(:,i) = (forwardfx-fx)/delta;
        end
end
end