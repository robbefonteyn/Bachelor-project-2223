function dist = dist_misval(A,B,MODE)

% dist = dist_misval(A,B,MODE)

nr_mis_val = sum(sum(~finite(A))) + sum(sum(~finite(B)));

dif=A-B*ones(1,size(A,2));
if nr_mis_val == 0
   if MODE==1
      dist=max(abs(dif),[],1);

   elseif MODE==2
      dist=sqrt(sum(dif.^2,1));

   else
      dist=sum(abs(dif),1);
   end
else
   if MODE==1
      dist=max(abs(dif),[],1);
   elseif MODE==2
      # Calculate the number of dimensions of the dif variable
      Nrdim=size(dif,1);
      # calculate the number of finite values (so the non missing values?)
      # look up what the 1 does after the comma
      Q=sum(finite(dif),1);
      # Set every place where the array Q has a value of 0 to NaN
      Q(find(Q==0))=NaN;
      # Set all infinite values in dif to 0
      dif(find(~finite(dif))) = 0;
      # Calculate the euclidian distance between (?)
      dist=sqrt((Nrdim*Q.^-1).*sum(dif.^2,1));
   else
      # Same calculation as the elseif part of this if statement,
      # move the following 4 lines under the else of the main if statement start
      Nrdim=size(dif,1);
      Q=sum(finite(dif),1);
      Q(find(Q==0))=NaN;
      dif(find(~finite(dif)))=0;
      # Only difference here in the calculation is that it's not the sqrt
      # and the sum of all absolute values is calculated in the second term
      dist=(Nrdim*Q.^-1).*sum(abs(dif),1);
   end
end
