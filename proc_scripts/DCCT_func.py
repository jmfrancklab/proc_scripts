from pylab import *
from pyspecdata import *
import matplotlib.lines as lines
from matplotlib.patches import FancyArrow, FancyArrowPatch
from pyspecdata.plot_funcs.image import imagehsv

label_spacing_multiplier = 50
allow_for_text_default = 10
allow_for_ticks_default = 70
text_height = 20

def DCCT(this_nddata,this_fig_obj,x=[],y=[],custom_scaling=False,**kwargs):
    this_nddata = this_nddata.C
    my_data = this_nddata.C
    my_data.human_units()
    print("DIMLABELS ARE",my_data.dimlabels)
    grid_bottom = 0.0
    bottom_pad = 0.15
    grid_bottom += bottom_pad
    grid_top = 1.0
    top_pad = 0.3
    grid_top -= top_pad
    total_spacing = 0.2
    a_shape = ndshape(this_nddata)
    num_dims = len(a_shape.dimlabels[:-2])
    divisions = []
    # should be looping in backward order from printed shape
    for j,thisdim in enumerate(a_shape.dimlabels[::-1][2:]):
        old = [j/2.0 for j in divisions]
        divisions = (old + [1])*(a_shape[thisdim]-1)+old
        print("for",thisdim,"I get",divisions)
    divisions = [j*total_spacing/sum(divisions) for j in divisions]
    axes_height = (grid_top-grid_bottom-total_spacing)/prod(a_shape.shape[:-2])
    axes_bottom = np.cumsum([axes_height+j for j in divisions]) # becomes ndarray
    axes_bottom = r_[0,axes_bottom]
    axes_bottom += grid_bottom
    axes_top = grid_bottom + grid_top
    fig = this_fig_obj
    ax_list = []
    yMajorLocator = lambda: mticker.MaxNLocator(steps=[1,2,5,10])
    majorLocator = lambda: mticker.MaxNLocator(min_n_ticks=4, steps=[1,2,5,10])
    minorLocator = lambda: mticker.AutoMinorLocator(n=5)

    LHS_pad = 0.01
    RHS_pad = 0.05
    LHS_labels = 0.13*num_dims
    width = 1.-(LHS_pad+RHS_pad+LHS_labels)

    for j,b in enumerate(axes_bottom):
        ax_list.append(axes([LHS_labels+LHS_pad,b,width,axes_height])) # lbwh
        if j == 0:
            print("OK")
            #AG: I ran the code with these 3 commented out and it still looks good
            #are these necessary or can we delete them?
            ax_list[-1].xaxis.set_major_locator(majorLocator())
            ax_list[-1].xaxis.set_minor_locator(minorLocator())
            ax_list[-1].set_ylabel(None)
        elif (j == len(axes_bottom)-1):
            ax_list[-1].xaxis.set_major_locator(majorLocator())
            ax_list[-1].set_xlabel(None)
            #for the minor ticks, use no labels; default NullFormatter
            ax_list[-1].xaxis.set_minor_locator(minorLocator())
            ax_list[-1].xaxis.tick_top()
            labels = [item.get_text() for item in ax_list[-1].get_xticklabels()]
            empty_string_labels = ['']*len(labels)
            ax_list[-1].set_xticklabels(empty_string_labels)
            ax_list[-1].set_xlabel(None)
        else:
            ax_list[-1].xaxis.set_ticks([])
            ax_list[-1].get_xaxis().set_visible(False)
            ax_list[-1].set_xlabel(None)
        ax_list[-1].set_ylabel(a_shape.dimlabels[-2])
        ax_list[-1].yaxis.set_minor_locator(minorLocator())
        ax_list[-1].yaxis.set_ticks_position('both')
        for tick in ax_list[-1].get_yticklabels():
            tick.set_rotation(0)

    if len(a_shape.dimlabels) > 3:
        A = this_nddata.smoosh(a_shape.dimlabels[:-2],'smooshed',noaxis=True)
        A.reorder('smooshed',first=True)
    else:
        A = this_nddata.C
        A_data = A.C
        A.rename(a_shape.dimlabels[:-2][0],'smooshed')

    def draw_span(ax1, ax2, label, this_label_num,
            allow_for_text=allow_for_text_default, allow_for_ticks=allow_for_ticks_default):
        x1,y1 = ax1.transAxes.transform(r_[0,1])
        x2,y2 = ax2.transAxes.transform(r_[0,0])
        x1-=allow_for_ticks
        x_text = x1-allow_for_text
        x2-=allow_for_ticks
        # following line to create an offset for different dimension labels
        label_spacing = this_label_num*label_spacing_multiplier
        x1,y1 = fig.transFigure.inverted().transform(r_[x1-label_spacing,y1])
        x_text,_ = fig.transFigure.inverted().transform(r_[x_text-label_spacing,0])
        x2,y2 = fig.transFigure.inverted().transform(r_[x2-label_spacing,y2])
        lineA = lines.Line2D([x1,x2],[y1,y2],
                linewidth=3, color='k', transform=fig.transFigure,
                clip_on=False)
        text(x_text, (y2+y1)/2, label, va='center', ha='right', rotation=90, transform=fig.transFigure, color='k')
        fig.add_artist(lineA)

    label_placed = zeros(num_dims)

    def place_labels(ax1, label, label_placed, this_label_num, check_for_label_num = True,
            allow_for_text=allow_for_text_default,
            allow_for_ticks=allow_for_ticks_default, y_space=30, arrow_width_px=3,
            push_bottom_label=0):
        if check_for_label_num:
            if not label_placed[this_label_num]:
                x1disp,y1disp = ax1.transAxes.transform(r_[0,0])
                label_spacing = this_label_num*label_spacing_multiplier
                x_text_disp = x1disp-(allow_for_text+allow_for_ticks)-label_spacing-text_height/2
                x_text,y1 = fig.transFigure.inverted().transform(r_[
                    x_text_disp+push_bottom_label,
                    y1disp-y_space])
                # push the arrow slightly below the topline of the text
                x_arrow,y_arrow = fig.transFigure.inverted().transform(r_[
                    x_text_disp,
                    y1disp-y_space-2*arrow_width_px])
                arrow_width,_ = fig.transFigure.inverted().transform(r_[arrow_width_px,0])
                dx,dy = fig.transFigure.inverted().transform(r_[
                    0,y_space-arrow_width_px])
                a = FancyArrow(x_arrow-arrow_width/2, y_arrow, dx, dy,
                        arrow_width,  alpha=0.1,  color='k')
                # could do fancier w/ the following, but need to mess w/ width parameters
                #arrow_base = r_[x_arrow-arrow_width/2, y_arrow]
                #a = FancyArrowPatch(arrow_base, arrow_base+r_[dx, dy],
                #        arrowstyle='|-|',
                #        alpha=0.1,  color='k')
                fig.add_artist(a)
                # the labels of the outer dimensions
                text(x_text, y1, label, va='top', ha='right', rotation=45,
                        transform=fig.transFigure, color='k')
                label_placed[this_label_num] = 1
        else:
            x1,y1disp = ax1.transAxes.transform(r_[0,0])
            label_spacing = this_label_num*80
            # from here https://stackoverflow.com/questions/44012436/python-matplotlib-get-position-of-xtick-labels
            # then searching for BBox docs
            print("tick locations",[j.get_window_extent().bounds for j in ax1.get_yticklabels()])
            x_textdisp = [j.get_window_extent().bounds for j in ax1.get_yticklabels()][0][0]
            x_textdisp -= text_height/2
            x_text,y1 = fig.transFigure.inverted().transform(r_[x_textdisp,y1disp-y_space])
            text(x_text, y1, my_data.unitify_axis(my_data.dimlabels[-2]), 
                    va='top', ha='right', rotation=45, transform=fig.transFigure, color='k')
            # push the arrow slightly below the topline of the text
            x_arrow,y_arrow = fig.transFigure.inverted().transform(r_[
                x_textdisp,
                y1disp-y_space-2*arrow_width_px])
            arrow_width,_ = fig.transFigure.inverted().transform(r_[arrow_width_px,0])
            dx,dy = fig.transFigure.inverted().transform(r_[
                0,y_space-arrow_width_px])
            a = FancyArrow(x_arrow-arrow_width/2,y_arrow,dx,dy,arrow_width, alpha=0.1, color='k')
            fig.add_artist(a)
            labels = [item.get_text() for item in ax1.get_xticklabels()]
            empty_string_labels = ['']*len(labels)
            ax1.set_xticklabels(empty_string_labels)
            ax1.set_xlabel(None)
            ax1.tick_params(bottom=False)
            ax1.set_ylabel(None)
            ax1.set_yticklabels(empty_string_labels)
            ax1.tick_params(left=False)

    imagehsvkwargs = {}
    for k,v in list(kwargs.items()):
        if k in ['black','logscale']:
            imagehsvkwargs[k] = kwargs.pop(k)
    

    for j in range(len(ax_list)):
        spacing,ax,x_first,origin,renumber = process_kwargs([('spacing',1),
            ('ax',ax_list[j]),
            ('x_first',False),
            ('origin','lower'),
            ('renumber',None)],kwargs,
            pass_through=True)
        if isinstance(x, list):
            x = np.array(my_data.getaxis('t2'))
        if isinstance(y, list):
            y = np.array(my_data.getaxis(my_data.dimlabels[-2]))
        if len(x)==0:
            x = [1,A.data.shape[1]]
        else:
            x = x.flatten()
        if len(y)==0:
            y = [1,A.data.shape[0]]
        else:
            y = y.flatten()
        dx = (x[-1]-x[0])/len(x)
        dy = (y[-1]-y[0])/len(y)
        if origin == 'lower':
            myext = (x[0]-dx/2.,x[-1]+dx/2.,y[0]-dy/2.,y[-1]+dy/2.)
        elif origin == 'upper':
            myext = (x[0]-dx/2.,x[-1]+dx/2.,y[-1]+dy/2.,y[0]-dy/2.)
        elif origin == 'flip':
            # {{{ need to flip
            myext = (x[-1]+dx/2.,x[0]-dx/2.,y[-1]+dy/2.,y[0]-dy/2.)
            # }}}
        else:
            raise ValueError("I don't understand the value you've set for the origin keyword argument")
        kwargs['origin'] = origin# required so that imshow now displays the image correctly
        
        if custom_scaling:
            scaling = 60.6856
            K = imagehsv(A['smooshed',j].data,**imagehsvkwargs,scaling=scaling)
        if not custom_scaling:
            K = imagehsv(A['smooshed',j].data,**imagehsvkwargs,scaling=abs(A).data.max())
        sca(ax_list[j])
        imshow(K,extent=myext,**kwargs)
        ax_list[j].set_ylabel(None)
        if not j == 0:
            ax_list[j].set_xlabel(None)
    # to drop into ax_list, just do
    # A.smoosh(a_shape.dimlabels, 'smooshed', noaxis=True)
    # in ax_list[0] put A['smooshed',0], etc
    idx = nddata(r_[0:prod(a_shape.shape[:-2])],[-1],['smooshed'])
    idx.chunk('smooshed',a_shape.dimlabels[:-2],a_shape.shape[:-2])
    remaining_dim = a_shape.dimlabels[:-2]
    depth = num_dims
    def decorate_axes(idx,remaining_dim,depth):
        thisdim=remaining_dim[0]
        print("This is remaining dim",remaining_dim)
        print("This dim is",thisdim)
        print(ndshape(idx))
        depth -= 1
        for j in range(a_shape[thisdim]):
            idx_slice = idx[thisdim,j]
            print("For",thisdim,"element",j,idx_slice.data.ravel())
            first_axes = ax_list[idx_slice.data.ravel()[0]]
            last_axes = ax_list[idx_slice.data.ravel()[-1]]
            draw_span(last_axes,first_axes,"%d"%(j),
                    this_label_num=depth)
            place_labels(ax_list[0],"%s"%my_data.unitify_axis('%s'%thisdim), label_placed,
                    this_label_num=depth)
            new_remaining_dim = remaining_dim[1:]
            if len(remaining_dim) > 1:
                decorate_axes(idx_slice,new_remaining_dim,depth)
    print("call recursive function")
    decorate_axes(idx,remaining_dim,depth)
    place_labels(axes([LHS_labels+LHS_pad,axes_bottom[0],width,0]),
            "%s"%(a_shape.dimlabels[-2]), label_placed,this_label_num=depth-1, 
            check_for_label_num = False, allow_for_text = -50)
    axes([LHS_labels+LHS_pad,axes_bottom[0],
        width,0]).set_xlabel(my_data.unitify_axis(my_data.dimlabels[-1]),
                labelpad=20)
    return 
