function mn = mean_misval(A,MODE)
# Used to calculate the mean of a given array and 2 options for the mode (1,2 or None)
# Can be used even if there are missing values (inf) in the array

# Calculate the number of missing values in the array A, same way as is done in dist_misval()
nr_mis_val = sum(sum(~finite(A)));

# If there are none, just calculate the mean of the data in the array depending on the mode
if nr_mis_val == 0
   if MODE==1
      # min(A,[],2) is a column vector with all the minimal values of each row on it
      # This just calculates the mean of the minimal value + the max value divided by 2
      # It is done to calculate the "mid-range" of the array
      mn = (min(A, [], 2) + max(A, [], 2)) / 2;

   elseif MODE == 2
      # Calculates the regular mean of the array
      mn = mean(A, 2);
   else
      # calculates the median as the mean of the array? when mode is set anything but 2 or 1
      mn=median(A,2);
   end
else

   if MODE==1
      # same exact calculation if nr_mis_val == 0 and mode == 1
      # (can be added to first if statement maybe, to avoid repitition?)
      # Calculates the mid-range again of the dataset
      mn=(min(A,[],2)+max(A,[],2))/2;

   elseif MODE==2
      # Q is a column vector with the sum of all the cells that have finite values for each row
      # (meaning that it counts the amount of 1's that isfinite() returns on each row)
      Q=sum(finite(A),2);
      # Replaces every 0 in Q (meaning that there are no finite values on that specific row of the array)
      Q(find(Q==0))=NaN;
      # Replaces the Inf values in the original array with the value 0
      A(find(~finite(A)))=0;
      # The dot operator here causes each value to be divided by the corresponding value in the Q column vector.
      # So mn is also a column vector as well
      mn = sum(A, 2)./Q;

   else
      # Nrrow is just the number of rows in the array
      Nrrow=size(A,1);
      for row=1:Nrrow
         # assign a median of a row to M calculated by leaving just the finite values of each row and using this
         # and if M is not empty (empty could be when a row is just Inf) it will add it to
         # mn on that same row. It adds NaN to mn if M isempty
         # The median is used as the mean when the data follows a normal distribution
         M=median(A(row,find(finite(A(row,:)))));
         if ~isempty(M)
            mn(row,1)=M;
         else
            mn(row,1)=NaN;
         end
      end
   end
end

