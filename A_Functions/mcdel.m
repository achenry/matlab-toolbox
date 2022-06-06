function [DEL, Cycles] = mcdel(time, data, m)
    
    %returns rainflow cycles and DELs, same results as MCrunch
    %this code is just a copy of the relevant parts of MCrunch code
    %F. Dunne 3/30/2011

    UCMult=0.5; %unclosed cycle multiplier (0 to 1)
    LUlt=25e6;   %ultimate load
    LMF=0;      %irrelevant when LUlt==inf
    
    %make sure data is in a vertical vector
    [a,b] = size(data);
    if a==1
        if b==1
            fprintf( '\n        Error: Data should be a vector.');
            Cycles  = [];
            DEL = [];
            return
        else
            data=data.';
        end
    else
        if b==1
        else
            fprintf( '\n        Error: Data should be a vector.');
            Cycles  = [];
            DEL = [];
            return
        end
    end
    %end make sure data is in a vertical vector
            
    
   % Finds rainflow cycles in the time series.
   %
   % Cycles(1) = Cycle range for variable means.
   % Cycles(2) = Cycle means.
   % Cycles(3) = Cycle range for fixed means.
   % Cycles(4) = Effective cycle weight.  Use one for full cycles, use UCMult for unclosed cycles.


      
         % Identify peaks and troughs in the time series.  The first and last points are considered peaks or troughs.
         % Sometimes the peaks can be flat for a while, so we have to deal with that nasty situation.
      [ Peaks, NumPeaks, PeakTimes ] = GetPeaks( data , time );


      % See if we have at least three points in the Peaks array.
      if ( NumPeaks < 3 )
         fprintf( '\n        This data has only %d peaks, so rainflow analysis is not possible.\n\n', NumPeaks );
         Cycles  = [];
         DEL = [];
         return
      end % if ( NumPeaks < 3 )


      % Identify the closed and unclosed cycles.  All cycle ranges, means, and weights are returned in the Cycles array.
      [Cycles] = GenCycles( Peaks, UCMult, LUlt, LMF );

      % Compute the simple, damage-equivalent loads.
      if(length(time)~=1)
          TimeFact=1/(time(end)-time(1));
      else
          TimeFact=1/time;
      end
      DEL = ( TimeFact*sum( Cycles(:,4).*( Cycles(:,3).^m ) ) )^( 1/m );
      Cycles=Cycles(:,3).';
     

      return
%=======================================================================
      function [Cycles] = GenCycles( Peaks, UCMult, LUlt, LMF)
      % Generate rainflow cycles.

      % Algorithm obtained from:

      %     Ariduru, Seçil (2004).  "Fatigue Life Calculation by Rainflow Cycle Counting Method."
      %     M.S. Thesis.  Ankara, Turkey: Middle East Technical University.

      % The example used in Section 3.2 of the thesis was used to debug this routine.
      % This routine also gives the exact same answers as Crunch.


            % Process the peaks and valleys.

         NumPeaks    = size( Peaks, 1 );
         RemainPeaks = zeros( NumPeaks, 1 );
         Ind         = 0;
         LenUC       = 0;
         NumCycles   = 0;
         Cycles      = zeros( int32( size( Peaks, 1 )/2 - 0.5 ), 4 );
         LFMargin    = LUlt - LMF;

         while ( true )


            if ( Ind < NumPeaks )

               Ind = Ind + 1;

               if ( Ind > NumPeaks ), break, end

               LenUC              = LenUC + 1;
               RemainPeaks(LenUC) = Peaks(Ind);

            end % if ( Ind < NumPeaks )


               % Make sure we have at least three peaks in the RemainPeaksing array.

            while ( LenUC < 3 )

               Ind = Ind + 1;

               if ( Ind > NumPeaks ), break, end

               LenUC              = LenUC + 1;
               RemainPeaks(LenUC) = Peaks(Ind);

            end % while ( LenUC < 3 )


               % Compute the newest and oldest active ranges.

            OldRange = abs( RemainPeaks(LenUC-1) - RemainPeaks(LenUC-2) );
            NewRange = abs( RemainPeaks(LenUC  ) - RemainPeaks(LenUC-1) );


               % If the new range is as large as the oldest active range, we found a cycle.
               % If LenUC is 3, it's a half cycle.  Add it to the list of cycles.

            while ( NewRange >= OldRange )

               NumCycles = NumCycles + 1;

               Cycles(NumCycles,1) = OldRange;
               Cycles(NumCycles,2) = 0.5*( RemainPeaks(LenUC-1) + RemainPeaks(LenUC-2) );
               if LUlt==inf
                   Cycles(NumCycles,3) = Cycles(NumCycles,1);
               else
                   Cycles(NumCycles,3) = Cycles(NumCycles,1)*LFMargin/( LUlt - abs( Cycles(NumCycles,2) ) );
               end
               if ( LenUC > 3 )
                  Cycles(NumCycles,4)              = 1.0;
                  RemainPeaks((LenUC-2):(LenUC-1)) = [];
                  LenUC                            = LenUC - 2;
               else
                  Cycles(NumCycles,4)    = UCMult;
                  RemainPeaks((LenUC-2)) = [];
                  LenUC                  = LenUC - 1;
               end % if ( LenUC > 3 )

               if ( LenUC >= 3 )
                  OldRange = abs( RemainPeaks(LenUC-1) - RemainPeaks(LenUC-2) );
                  NewRange = abs( RemainPeaks(LenUC  ) - RemainPeaks(LenUC-1) );
               else
                  NewRange = -1;
               end % if ( LenUC >= 3 )
               

            end % while ( NewRange >= OldRange )

            if ( Ind == NumPeaks ), break, end

         end % while


            % Add the unclosed cycles to the end of the Cycles matrix if the weight is not zero.

         if ( ( LenUC > 1 ) && ( UCMult > 0 ) )

            for Cyc=1:LenUC-1
               Cycles(NumCycles+Cyc,1) = abs ( RemainPeaks(Cyc) - RemainPeaks(Cyc+1) );
               Cycles(NumCycles+Cyc,2) = 0.5*( RemainPeaks(Cyc) + RemainPeaks(Cyc+1) );
               if LUlt==inf
                   Cycles(NumCycles+Cyc,3) = Cycles(NumCycles+Cyc,1);
               else
                   Cycles(NumCycles+Cyc,3) = Cycles(NumCycles+Cyc,1)*LFMargin/( LUlt - abs( Cycles(NumCycles+Cyc,2) ) );
               end
               Cycles(NumCycles+Cyc,4) = UCMult;
            end % for Cyc

         else

            LenUC = 1;

         end % if ( ( LenUC > 1 ) && ( UCMult > 0 ) )


            % Truncate the unused portion of the array.

         TotCycles = NumCycles + LenUC - 1;
         Cycles    = Cycles(1:TotCycles,:);


      end % function Cycles = GenCycles( Peaks, UCMult )
%=======================================================================
      function [ Peaks, NumPeaks , PeakTimes ] = GetPeaks( TimeSer , time )
      % Identify peaks and troughs in the time series.  The first and last points are considered peaks or troughs.
      % Sometimes the peaks can be flat for a while, so we have to deal with that nasty situation.

         Peaks    = TimeSer;
         NumPeaks = 1;
         LastDiff = 1;
         TSlen    = size( TimeSer, 1 );
         PeakTimes= time;

         for Pt=2:(TSlen-1)

            if ( TimeSer(Pt) == TimeSer(Pt+1) )                                                               % Is slope zero?  Don't update LastDiff is so.
               continue;
            elseif ( ( sign( TimeSer(Pt) - TimeSer(LastDiff) ) + sign( TimeSer(Pt+1) - TimeSer(Pt) ) ) ==  0 )    % Did slope change sign?
               NumPeaks        = NumPeaks + 1;
               Peaks(NumPeaks) = TimeSer(Pt);
               PeakTimes(NumPeaks) = time(Pt);
            end % if

            LastDiff = Pt;

         end % for Pt


            % Add the last point of the time series to the list of peaks.

         Peaks    = [ Peaks(1:NumPeaks); TimeSer(TSlen) ];
         PeakTimes= [ PeakTimes(1:NumPeaks);time(TSlen) ];
         NumPeaks = NumPeaks + 1;

      end % function [ Peaks, NumPeaks ] = GetPeaks( TimeSer )
end %mcdel

