close all
clear all
[num,txt,raw] = xlsread('100 Hz 1000 Samples.xlsx');

features_raw=num(1:end,7:1006);  
age=num(1:end,3:3);
fig, ax = plt.subplots(1, 1) 
ax.hist(age) 
   
  
%Set title 
ax.set_title("Title") 
  
% adding labels 
ax.set_xlabel('x-label') 
ax.set_ylabel('y-label') 
  
% Make some labels. 
rects = ax.patches 
labels = ["label%d" % i for i in range(len(rects))] 
  
%for rect, label in zip(rects, labels): 
    height = rect.get_height() 
    ax.text(rect.get_x() + rect.get_width() / 2, height+0.01, label, 
            ha='center', va='bottom') 
  
%Show plot 
plt.show() 