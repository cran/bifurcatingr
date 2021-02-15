#' Plotting Bifurcating Autoregressive Trees
#'
#' This function graphs bifurcating autoregressive (BFA) tree data.
#'
#' @param z a numeric vector containing the tree data
#' @param digits an integer indicating the number of decimal places to be
#'   displayed in vertex labels
#' @param shape the shape of the vertex. Currently “circle”, “square”,
#'   “csquare”, “rectangle”, “crectangle”, “vrectangle”, and “none” are
#'   supported. Defaults to “none” which does not display the vertices at all.
#' @param vertex.size a numeric scalar or vector defining the size of the vertex
#'   or vertices. If a vector is supplied, vertex sizes may differ. Defaults to
#'   10.
#' @param text.size the font size of vertex labels. Defaults to 1.
#' @param text.color the color of vertex labels. If it is a character vector,
#'   then it may either contain integer values, named colors or RGB specified
#'   colors with three or four bytes. Defaults to "black".
#' @param vertex.color the fill color of the vertex. If you don't want some or
#'   all vertices to have any color, supply NA. The default is "gold". See also
#'   the options in \code{text.color}.
#' @param arrow.size the size of the arrows. The default value is 0.5.
#' @param arrow.width the width of the arrows. The default value is 0.5.
#' @param arrow.color the color of the arrows. The default is "black". See also
#'   the options in \code{text.color}.
#' @param plot.margin the amount of empty space around the plot, it is a numeric
#'   vector of length four. Usually values between 0 and 0.5 are meaningful, but
#'   negative values are also possible and in that case it will make the plot
#'   zoom in to a part of the graph. If it is shorter than four, recycling will
#'   occur. The default value is -0.3.
#' @return A binary tree displaying the BFA data. 
#' @export
#' @details For more details about the graph options see
#'   \code{\link[igraph]{igraph.plotting}}.
#' @examples
#' z <- bfa.tree.gen(31, 1, 1, 1, 0.5, 0.5, 0, 10, c(0.7))
#' bfa.tree.plot(z)
#' bfa.tree.plot(z,shape= "circle")
#' bfa.tree.plot(z,shape= "circle", text.color="white", vertex.color = "darkgrey")
bfa.tree.plot<- function(z, digits, shape= "none", vertex.size = 10, text.size = 1,
                         text.color="black", vertex.color = "gold", arrow.size = 0.5,
                         arrow.width = 0.5, arrow.color="black", plot.margin = -0.3){
  tree <- igraph::make_tree(n=length(z), children=2)
  igraph::V(tree)$name <- round(z, digits=digits)
  lay <- igraph::layout_as_tree(tree,flip.y = TRUE)
  plot(tree, layout=lay, vertex.shape= shape, vertex.size=vertex.size, vertex.label.cex=text.size,
       vertex.label.color=text.color,vertex.color=vertex.color, edge.color=arrow.color,
       edge.arrow.size = arrow.size, edge.width = arrow.width, margin=plot.margin)
}

