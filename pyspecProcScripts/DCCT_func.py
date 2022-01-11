from pylab import *
from pyspecdata import *
import matplotlib.lines as lines
from matplotlib.patches import FancyArrow, FancyArrowPatch, Circle
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation,IdentityTransform
from pyspecdata.plot_funcs.image import imagehsv
import logging
@FuncFormatter
def ph2(x,pos):
    ordered_labels = {}
    n_ph = 2 
    this_max_coh_jump = 1
    all_possibilities = empty((int((2*this_max_coh_jump+1)/n_ph)+1)*n_ph)
    all_possibilities[:] = nan
    all_possibilities[:this_max_coh_jump+1] = r_[0:this_max_coh_jump+1]
    all_possibilities[-this_max_coh_jump:] = r_[-this_max_coh_jump:0]
    all_possibilities = all_possibilities.reshape((-1,n_ph))
    labels_in_order = []
    for j in range(n_ph):
        temp = all_possibilities[:,j]
        if j == 0:
            temp = ', '.join(['%d'%j for j in 
                temp[isfinite(temp)]])
        else:    
            temp = ', '.join(['%+d'%j for j in
                temp[isfinite(temp)]])
        if len(temp) == 0: temp = 'X'
        labels_in_order.append(temp)
    ordered_labels['ph2'] = labels_in_order
    if x == 0:
        temp = ("%s")%ordered_labels['ph2'][0]
    #else:
    #    temp = ('%s')%ordered_labels['ph2'][]
    return temp


diagnostic = False
def DCCT(this_nddata, this_fig_obj, x=[], y=[], custom_scaling=False,
        grid_bottom = 0.0,
        bottom_pad = 0.15,
        grid_top = 1.0,
        top_pad = 0.1,
        total_spacing = 0.055,
        label_spacing_multiplier = 50,
        allow_for_text_default = 10,
        allow_for_ticks_default = 70,
        text_height = 40,
        LHS_pad = 0.01,
        RHS_pad = 0.05,
        shareaxis = False,
        cmap = None,
        pass_frq_slice = False,
        frq_slice = [],
        just_2D = False,
        scaling_factor=1,
        max_coh_jump={'ph1':2,'ph2':1},
        direct='t2',
        plot_title = 'DCCT',
        **kwargs):
    """DCCT plot

    Parameters
    ==========
    shareaxis: boolean
        subplots scale together, but currently, this means there must be tick labels on both top and bottom
    """
    my_data = this_nddata.C
    ordered_labels = {}
    # {{{ Generate alias labels - goes to scientific fn
    for this_dim in [j for j in my_data.dimlabels if j.startswith('ph')]:
        n_ph = ndshape(my_data)[this_dim]
        this_max_coh_jump = max_coh_jump[this_dim]
        all_possibilities = empty((int((2*this_max_coh_jump+1)/n_ph)+1)*n_ph)
        all_possibilities[:] = nan
        all_possibilities[:this_max_coh_jump+1] = r_[0:this_max_coh_jump+1]
        all_possibilities[-this_max_coh_jump:] = r_[-this_max_coh_jump:0]
        all_possibilities = all_possibilities.reshape((-1,n_ph))
        labels_in_order = []
        for j in range(n_ph):
            temp = all_possibilities[:,j]
            if j == 0:
                temp = ', '.join(['%d'%j for j in 
                    temp[isfinite(temp)]])
            else:    
                temp = ', '.join(['%+d'%j for j in
                    temp[isfinite(temp)]])
            if len(temp) == 0: temp = 'X'
            labels_in_order.append(temp)
        ordered_labels[this_dim] = labels_in_order
    # }}}
    real_data = False
    if cmap is not None:
        assert all(isclose(my_data.data.imag,0)), "In order to use a color map, you must pass real data"
        if type(cmap) == str:
            cmap = get_cmap(cmap)
            my_data.data = my_data.data.real
            real_data = True
    my_data.human_units()
    grid_bottom += bottom_pad
    grid_top -= top_pad
    a_shape = ndshape(this_nddata)
    num_dims = len(a_shape.dimlabels[:-2])
    divisions = []
    
    # should be looping in backward order from printed shape
    for j,thisdim in enumerate(a_shape.dimlabels[::-1][2:]):
        old = [j/2.0 for j in divisions]
        divisions = (old + [1])*(a_shape[thisdim]-1)+old
        logging.debug(strm("for",thisdim,"I get",divisions))
    divisions = [j*total_spacing/sum(divisions) for j in divisions]
    axes_height = (grid_top-grid_bottom-total_spacing)/prod(a_shape.shape[:-2])
    axes_bottom = np.cumsum([axes_height+j for j in divisions]) # becomes ndarray
    axes_bottom = r_[0,axes_bottom]
    axes_bottom += grid_bottom
    axes_top = grid_bottom + grid_top
    fig = this_fig_obj
    ax_list = []
    yMajorLocator = lambda: mticker.MaxNLocator(nbins='auto',
                        steps=[1,2,5,10])
    majorLocator = lambda: mticker.MaxNLocator(nbins='auto',
            steps=[1,2,2.5,5,10])
    minorLocator = lambda: mticker.AutoMinorLocator(n=5)
    LHS_labels,_ = fig.transFigure.inverted().transform(
            (label_spacing_multiplier*num_dims + allow_for_ticks_default, 0))
    width = 1.-(LHS_pad+RHS_pad+LHS_labels)
    for j,b in enumerate(axes_bottom):
        if j !=0 and shareaxis:
            ax_list.append(axes([LHS_labels+LHS_pad,b,width,axes_height],
                sharex=ax_list[0],
                sharey=ax_list[0],
                )) # lbwh
        else:
            ax_list.append(axes([LHS_labels+LHS_pad,b,width,axes_height])) # lbwh
    # {{{ adjust tick settings -- AFTER extents are set
    # {{{ bottom subplot
    ax_list[0].xaxis.set_major_locator(majorLocator())
    ax_list[0].xaxis.set_minor_locator(minorLocator())
    ax_list[0].set_ylabel(None)
    ax_list[0].set_xlabel(my_data.unitify_axis(my_data.dimlabels[-1]),
                labelpad=20)
    # }}}
    # {{{ intermediate subplots
    for j in range(1,len(axes_bottom)-1):
        ax_list[j].xaxis.set_ticks([])
        ax_list[j].get_xaxis().set_visible(False)
        ax_list[j].set_xlabel(None)
    # }}}
    # {{{ top subplot
    ax_list[-1].xaxis.set_major_locator(majorLocator())
    #for the minor ticks, use no labels; default NullFormatter
    ax_list[-1].xaxis.set_minor_locator(minorLocator())
    ax_list[-1].xaxis.tick_top()
    if not shareaxis:
        ax_list[-1].set_xticklabels([])
    # }}}
    # {{{ all subplots
    for j in range(0,len(axes_bottom)):
        ax_list[j].set_ylabel(a_shape.dimlabels[-2])
        inner_dim = a_shape.dimlabels[-2]
        inner_dim = str(inner_dim)
        if inner_dim == 'ph2':
            logging.debug("Inner dimension is phase cycling dimension")
            ax_list[j].yaxis.set_major_formatter(ph2)
            ax_list[j].yaxis.set_major_locator(MaxNLocator(integer=True))
        else:
            ax_list[j].yaxis.set_minor_locator(minorLocator())
            ax_list[j].yaxis.set_ticks_position('both')
        for tick in ax_list[j].get_yticklabels():
            tick.set_rotation(90)
            tick.set_va('center')
    # }}}
    # }}}

    if len(a_shape.dimlabels) > 3:
        A = this_nddata.C.smoosh(a_shape.dimlabels[:-2],'smooshed',noaxis=True)
        A.reorder('smooshed',first=True)
    else:
        A = this_nddata.C
        A_data = A.C
        A.rename(a_shape.dimlabels[:-2][0],'smooshed')

    def draw_span(ax1, ax2, label, this_label_num,
            allow_for_text=allow_for_text_default, 
            allow_for_ticks=allow_for_ticks_default):
        x1,y1 = ax1.transAxes.transform(r_[0,0.95])
        x2,y2 = ax2.transAxes.transform(r_[0,0.05])
        x1-=allow_for_ticks
        x_text = x1-allow_for_text
        x2-=allow_for_ticks
        # following line to create an offset for different dimension labels
        label_spacing = this_label_num*label_spacing_multiplier
        x1,y1 = fig.transFigure.inverted().transform(r_[x1-label_spacing,y1])
        x_text,_ = fig.transFigure.inverted().transform(r_[x_text-label_spacing,0])
        x2,y2 = fig.transFigure.inverted().transform(r_[x2-label_spacing,y2])
        lineA = lines.Line2D([x1,x2],[y1,y2],
                linewidth=1, color='k', transform=fig.transFigure,
                clip_on=False)
        text(x_text, (y2+y1)/2, label, va='center', ha='right', rotation=90, transform=fig.transFigure, color='k')
        fig.add_artist(lineA)

    label_placed = zeros(num_dims)

    def place_labels(ax1, label, label_placed, this_label_num, check_for_label_num = True,
            allow_for_text=allow_for_text_default,
            allow_for_ticks=allow_for_ticks_default,
            y_space_px=30, arrow_width_px=4.0,
            arrow_head_vs_width=3):
        if not check_for_label_num or not label_placed[this_label_num]:
            x_axorigindisp,y_axorigindisp = ax1.transAxes.transform(r_[0,0])
            # {{{ determine the x and y position of the label in display coords
            if check_for_label_num:
                # the labels of the outer dimensions
                label_spacing = this_label_num*label_spacing_multiplier
                Dx_textdisp = -(allow_for_text+allow_for_ticks)-label_spacing-text_height/2
                x_textdisp = x_axorigindisp+Dx_textdisp
            else:
                # same as above, but determine text
                # position based on tick labels
                label = my_data.unitify_axis(my_data.dimlabels[-2])
                x_axorigindisp,y_axorigindisp = ax1.transAxes.transform(r_[0,0])
                # from here https://stackoverflow.com/questions/44012436/python-matplotlib-get-position-of-xtick-labels
                # then searching for BBox docs
                logging.debug(strm( "tick locations",[j.get_window_extent().bounds for j in ax1.get_yticklabels()] ))
                x_textdisp = [j.get_window_extent().bounds for j in ax1.get_yticklabels()][0][0]
            x_textdisp -= 3.7*text_height - arrow_width_px
            y_textdisp = -25.0 
            if diagnostic:
                a = Circle((x_textdisp, y_axorigindisp-y_textdisp), 3,
                        clip_on=False,
                        transform=None,
                        color='b')
                fig.add_artist(a)
            # }}}
            axis_to_figure = ax1.transAxes + fig.transFigure.inverted()
            ax_x,ax_y = axis_to_figure.transform(r_[0,0])
            AnArrow = FancyArrow(x_textdisp,y_textdisp,
                    2, 5,
                    width=arrow_width_px,
                    clip_on=False,
                    transform=(IdentityTransform()
                        +ScaledTranslation(ax_x,ax_y,fig.transFigure)),#fig.transFigure,
                    alpha=0.1,  color='k')
            # could do fancier w/ the following, but need to mess w/ width parameters
            #arrow_base = r_[x_arrowbase_fig-arrow_width/2, y_arrowbase_fig]
            #a = FancyArrowPatch(arrow_base, arrow_base+r_[dx, dy],
            #        arrowstyle='|-|',
            #        alpha=0.1,  color='k')
            fig.add_artist(AnArrow)
            if diagnostic:
                a = Line2D([x_arrowbase_fig, x_arrowbase_fig+dx],
                        [y_arrowbase_fig, y_arrowbase_fig+dy],
                        transform=fig.transFigure,
                        color='r')
                fig.add_artist(a)
                a = Circle((x_textfig, y_textfig), 0.01,
                        clip_on=False,
                        transform=fig.transFigure,
                        color='r')
                fig.add_artist(a)
            x_textfig = x_textdisp + arrow_width_px
            y_textfig = y_textdisp - 5.0
            text(x_textfig, y_textfig, label,
                    va='top', ha='right', rotation=45,
                    clip_on=False,
                    transform=(IdentityTransform()
                        +ScaledTranslation(ax_x,ax_y,fig.transFigure)),
                    color='k')
            if check_for_label_num:
                label_placed[this_label_num] = 1
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
            x = np.array(my_data.getaxis(direct))
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
        
        if real_data:
            kwargs['cmap'] = cmap
            K = A['smooshed',j].data / abs(A).data.max()
        else:
            if custom_scaling:
                #scaling = 60.6856
                K = imagehsv(A['smooshed',j].data,**imagehsvkwargs,scaling=scaling_factor)
            if not custom_scaling:
                K = imagehsv(A['smooshed',j].data,**imagehsvkwargs,scaling=abs(A).data.max())
        sca(ax_list[j])
        imshow(K,extent=myext,**kwargs)
        ax_list[j].set_ylabel(None)
        if pass_frq_slice:
            start_y = A.getaxis(A.dimlabels[1])[0]
            stop_y = A.getaxis(A.dimlabels[1])[-1]
            ax_list[j].fill([x[0],frq_slice[0],frq_slice[0],x[0]],
                [start_y,start_y,stop_y,stop_y],
                fill=None,alpha=0.7,hatch='//')
            ax_list[j].fill([frq_slice[-1],x[-1],x[-1],frq_slice[-1]],
                    [start_y,start_y,stop_y,stop_y],
                    fill=None,alpha=0.7,hatch='//')
    # to drop into ax_list, just do
    # A.smoosh(a_shape.dimlabels, 'smooshed', noaxis=True)
    # in ax_list[0] put A['smooshed',0], etc
    idx = nddata(r_[0:prod(a_shape.shape[:-2])],[-1],['smooshed'])
    idx.chunk('smooshed',a_shape.dimlabels[:-2],a_shape.shape[:-2])
    remaining_dim = a_shape.dimlabels[:-2]
    depth = num_dims
    def decorate_axes(idx,remaining_dim,depth):
        thisdim=remaining_dim[0]
        logging.debug(strm("This is remaining dim",remaining_dim))
        logging.debug(strm("This dim is",thisdim))
        depth -= 1
        for j in range(a_shape[thisdim]):
            idx_slice = idx[thisdim,j]
            logging.debug(strm("For",thisdim,"element",j,idx_slice.data.ravel()))
            first_axes = ax_list[idx_slice.data.ravel()[0]]
            last_axes = ax_list[idx_slice.data.ravel()[-1]]
            if my_data.get_ft_prop(thisdim) == True:    
                if j == 0:
                    draw_span(last_axes,first_axes,("%s")%ordered_labels[thisdim][0],this_label_num=depth)
                else:
                    draw_span(last_axes,first_axes,('%s')%ordered_labels[thisdim][j],
                        this_label_num=depth)
            else:
                draw_span(last_axes,first_axes,("%d")%(j),this_label_num=depth)

            place_labels(ax_list[0],"%s"%my_data.unitify_axis('%s'%thisdim), label_placed,
                    this_label_num=depth)
            new_remaining_dim = remaining_dim[1:]
            if len(remaining_dim) > 1:
                decorate_axes(idx_slice,new_remaining_dim,depth)
    decorate_axes(idx,remaining_dim,depth)
    place_labels(ax_list[0],
            "%s"%(a_shape.dimlabels[-2]), label_placed,this_label_num=depth-1, 
            check_for_label_num = False, allow_for_text = -50)
    plt.title(plot_title)
    if just_2D:
        return LHS_pad+LHS_labels,axes_bottom[0],width,axes_bottom[-1]-top_pad
    else:
        return ax_list, LHS_pad+LHS_labels,axes_bottom[-1]+axes_height,width,top_pad-RHS_pad

