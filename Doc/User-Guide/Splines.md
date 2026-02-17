# Splines {#splines}

Description of Splines and how to read them in.

MaCh3 supports both event-by-event and binned splines:

# Binned splines
TODO!!



# Event-by-event splines
Format to initlaise event-by-event `std::vector< std::vector<TSpline3_red*> > MasterSpline` with size `[NEvents][NSplineParams]`. If the given event doesn't have a spline then it should be set to `NULL`.

### SplineMonolith Fast Load File
To speed up the process of loading files there is a function `PrepareSplineFile()` which will dump spline information into a prepossed file with only relevant information. Then there is the prepared constructor `SMonolith::SMonolith(std::string FileName)`.

### GPU
Since the calcuatalion of event by event spline is quite often the bottleneck there is GPU support for these calculations. From developer and user point of view code works the same. If GPU support is enabled MaCh3 will copy spline information to GPU and perform calcaution there.

![image](https://github.com/mach3-software/MaCh3/assets/45295406/34f7aeaf-d40e-4e8c-96f7-a1eb197603ea)


### TF1
There is the possibility to use TF1 rather than cubic interpolation. Which helps significantly reduce RAM. However, it may not work for every spline. Therefore it is advised to check the LLH scan before making an analysis choice.
![image](https://github.com/mach3-software/MaCh3/assets/45295406/bf87261d-60cd-4a87-b83a-46c6afa3f89e)


# Available spline interpolation
At MaCh3 following interpolation are available for cubic splines:
<ol>
<li> TSpline3 </li>
<li> Linear </li>
<li> Monotonic </li>
<li> Akima </li>

<img width="1287" height="937" alt="image" src="https://github.com/user-attachments/assets/32fc41c0-3063-47c7-9353-fb24087213dc" />
<img width="1275" height="924" alt="image" src="https://github.com/user-attachments/assets/1edb571d-9f08-4fe7-8de2-4cac30a1fe95" />
<img width="1275" height="896" alt="image" src="https://github.com/user-attachments/assets/e6682bb2-193f-4cb8-8d65-b1780af9590d" />
