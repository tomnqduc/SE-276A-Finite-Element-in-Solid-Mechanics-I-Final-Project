function [nodesCoord, elemDict] = uniformRectangularMesh4Node(Lx, Ly, noElemX, noElemY)
    totalElem = noElemX * noElemY;
    totalNode = (noElemX+1) * (noElemY);
    nodePerElem = 4;
    x = linspace(0,Lx,noElemX+1);
    y = linspace(0,Ly,noElemY+1);
    x_mesh = meshgrid(x,y);
    y_mesh = meshgrid(y,x)';

    nodesCoord = zeros(totalNode,2); % Define coordinates of each node
    elemDict = zeros(totalElem,nodePerElem); % Define which node form an element
    
    currNode = 1;
    
    
    for y=1:noElemY+1
        for x=1:noElemX+1
            nodesCoord(currNode,1) = x_mesh(y,x);
            nodesCoord(currNode,2) = y_mesh(y,x);
            currNode = currNode + 1;
        end
    end
    
    for x=1:noElemX
        for y=1:noElemY
            currElemNumber = (y-1)*noElemX + x;
            if x==1 && y==1 %initiate the node of the 1st element
                elemDict(currElemNumber,1) = 1;
                elemDict(currElemNumber,2) = 2;
                elemDict(currElemNumber,3) = noElemX + 3;
                elemDict(currElemNumber,4) = noElemX + 2;

            elseif x==1
                elemDict(currElemNumber,1) = elemDict(currElemNumber-noElemX,4);
                elemDict(currElemNumber,2) = elemDict(currElemNumber-noElemX,3);
                elemDict(currElemNumber,3) = elemDict(currElemNumber,2) + noElemX+1;
                elemDict(currElemNumber,4) = elemDict(currElemNumber,1) + noElemX+1;
            else
                elemDict(currElemNumber,1) = elemDict(currElemNumber-1,2);
                elemDict(currElemNumber,4) = elemDict(currElemNumber-1, 3);
                elemDict(currElemNumber,3) = elemDict(currElemNumber,4) + 1;
                elemDict(currElemNumber,2) = elemDict(currElemNumber,1) + 1;
                
            end
        end
    end


end