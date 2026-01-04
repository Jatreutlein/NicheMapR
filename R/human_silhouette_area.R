#' human_silhouette_area
#'
#' Calculates whole-body silhouette area as a function of solar zenith and azimuth
#' angles using an endotherm geometric model, and compares it to the empirical
#' Underwood & Ward formulation. A polynomial correction is fitted to the ratio
#' between the two methods as a function of zenith angle.
#'
#' The body is decomposed into head, trunk, arms, and legs, with geometry
#' determined by mass, density, fat fraction, and shape parameters. Silhouette
#' area is summed across parts and optionally plotted.
#'
#' @encoding UTF-8
#'
#' @param ZENs numeric vector. Solar zenith angles (degrees).
#' @param AZIMUTHs numeric vector. Solar azimuth angles (degrees), same length as
#'   \code{ZENs}.
#' @param MASSs numeric vector (length 4). Masses of body parts (kg) in the order:
#'   head, trunk, arms, legs.
#' @param DENSITYs numeric vector (length 4). Tissue densities of body parts
#'   (kg m\eqn{^{-3}}).
#' @param INSDEPDs numeric vector (length 4). Insulation depth parameters (m) for
#'   each body part.
#' @param SHAPE_Bs numeric vector (length 4). Shape coefficients controlling body
#'   part geometry.
#' @param SUBQFATs numeric vector (length 4). Subcutaneous fat distribution
#'   factors (unitless).
#' @param FATPCTs numeric vector (length 4). Percent fat content of each body
#'   part.
#' @param PJOINs numeric vector (length 4). Fractional joint overlap parameters
#'   for body parts.
#' @param plot.sil logical. If \code{TRUE}, plots silhouette area comparisons and
#'   the fitted correction factor.
#'
#' @details
#' Silhouette area is computed using the \code{GEOM_ENDO} routine for each body
#' part and summed to obtain total projected area. The Underwood & Ward silhouette
#' approximation is evaluated for comparison, and their ratio is modelled as a
#' polynomial function of zenith angle for use in diffuse radiation corrections.
#'
#' @return A list with components:
#' \describe{
#'   \item{fit}{Linear model object fitting the Underwood/endoR silhouette ratio
#'   as a polynomial function of zenith angle.}
#'   \item{sil.underwood}{Numeric vector of silhouette areas (m\eqn{^{2}})
#'   calculated using the Underwood & Ward formulation.}
#'   \item{sil.HomoTherm}{Numeric vector of silhouette areas (m\eqn{^{2}})
#'   calculated using the endotherm geometric model.}
#'   \item{AREA}{Maximum total body surface area used for solar exchange
#'   (m\eqn{^{2}}).}
#' }
#'
#' @export
human_silhouette_area <- function(ZENs = seq(0, 90, 1), # degrees, zenith angles
                           AZIMUTHs = rep(0, length(ZENs)),
                           MASSs = c(5.32, 35.07, 3.43, 11.34), # kg, masses per part
                           DENSITYs = rep(1050, 4), # kg/m3, densities per part
                           INSDEPDs = c(1e-02, rep(6e-03, 3)),
                           SHAPE_Bs = c(1.6, 1.85, 12.5, 6.5),
                           SUBQFATs = rep(1, 4),
                           FATPCTs = c(5 * 0.1, 36 * 0.2, 10, 23),
                           PJOINs = c(0.02628205, 0.08049954, 0.01923077, 0.03333333),
                           plot.sil = TRUE){

  SHAPEs <- c(4, 1, 1, 1)
  ORIENTs <- rep(3, 4)
  for(j in 1:length(ZENs)){
    GEOM.head <- c(ZENs[j], GEOM_ENDO(MASSs[1], DENSITYs[1], DENSITYs[1], FATPCTs[1], SHAPEs[1], INSDEPDs[1], SUBQFATs[1], SHAPE_Bs[1], SHAPE_Bs[1], 0, 0, PJOINs[1], 0, ORIENTs[1], ZEN = ZENs[j]))
    GEOM.trunk <- c(ZENs[j], GEOM_ENDO(MASSs[2], DENSITYs[2], DENSITYs[2], FATPCTs[2], SHAPEs[2], INSDEPDs[2], SUBQFATs[2], SHAPE_Bs[2], SHAPE_Bs[2], 0, 0, PJOINs[2], 0, ORIENTs[2], ZEN = ZENs[j]))
    GEOM.arm <- c(ZENs[j], GEOM_ENDO(MASSs[3], DENSITYs[3], DENSITYs[3], FATPCTs[3], SHAPEs[3], INSDEPDs[3], SUBQFATs[3], SHAPE_Bs[3], SHAPE_Bs[3], 0, 0, PJOINs[3], 0, ORIENTs[3], ZEN = ZENs[j]))
    GEOM.leg <- c(ZENs[j], GEOM_ENDO(MASSs[4], DENSITYs[4], DENSITYs[4], FATPCTs[4], SHAPEs[4], INSDEPDs[4], SUBQFATs[4], SHAPE_Bs[4], SHAPE_Bs[4], 0, 0, PJOINs[4], 0, ORIENTs[4], ZEN = ZENs[j]))
    if(j == 1){
      GEOM.heads <- GEOM.head
      GEOM.trunks <- GEOM.trunk
      GEOM.arms <- GEOM.arm
      GEOM.legs <- GEOM.leg
    }else{
      GEOM.heads <- rbind(GEOM.heads, GEOM.head)
      GEOM.trunks <- rbind(GEOM.trunks, GEOM.trunk)
      GEOM.arms <- rbind(GEOM.arms, GEOM.arm)
      GEOM.legs <- rbind(GEOM.legs, GEOM.leg)
    }
  }
  GEOM.heads <- as.data.frame(GEOM.heads)
  GEOM.trunks <- as.data.frame(GEOM.trunks)
  GEOM.arms <- as.data.frame(GEOM.arms)
  GEOM.legs <- as.data.frame(GEOM.legs)

  GEOM.lab <- c("ZEN", "VOL", "D", "MASFAT", "VOLFAT", "ALENTH", "AWIDTH", "AHEIT", "ATOT", "ASIL", "ASILN", "ASILP", "GMASS", "AREASKIN", "FLSHVL", "FATTHK", "ASEMAJ", "BSEMIN", "CSEMIN", "CONVSK", "CONVAR", "R1", "R2")

  colnames(GEOM.heads) <- GEOM.lab
  colnames(GEOM.trunks) <- GEOM.lab
  colnames(GEOM.arms) <- GEOM.lab
  colnames(GEOM.legs) <- GEOM.lab
  AREA <- max(GEOM.heads$ATOT + GEOM.trunks$ATOT + GEOM.arms$ATOT * 2 + GEOM.legs$ATOT * 2)

  # compare with Underwood and Ward calculation
  thetas <- (90 - ZENs) * pi / 180
  psis <- AZIMUTHs * pi / 180
  sil.underwood <- (0.043 * sin(thetas) + 2.997 * cos(thetas) * (0.02133 * cos(psis) ^ 2 + 0.0091 * sin(psis)^2) ^ 0.5) / 1.81 * AREA
  sil.endoR <- GEOM.heads$ASIL + GEOM.trunks$ASIL + GEOM.arms$ASIL * 2 + GEOM.legs$ASIL * 2
  if(plot.sil){
    par(mfrow = c(2, 1))
    par(oma = c(2, 1, 2, 2) + 0.1)
    par(mar = c(3, 3, 1.5, 1) + 0.1)
    par(mgp = c(2, 1, 0))
    plot(ZENs, sil.endoR, type = 'l', ylim = c(0, max(sil.endoR)), ylab = 'silhouette area, m2', xlab = 'zenith angle, degrees')
    points(ZENs, sil.underwood, type = 'l', lty = 2)
    legend(c(0, max(sil.endoR)), cex = 0.75, legend = c('endoR', 'Underwood & Wang'), lty = c(1, 2), bty = 'n')
  }
  # model difference as a function of zenith and adjust diffuse solar fraction
  delta.calc <- sil.underwood / sil.endoR
  fit <- lm(delta.calc ~ poly(ZENs, 8))
  if(plot.sil){
    plot(ZENs, delta.calc, type = 'l', ylab = 'Underwood/endoR silhouette area', xlab = 'zenith angle, degrees')
    points(predict(fit), type = 'l', lty = 2)
    legend(max(ZENs) * 0.5, .7, cex = 0.75, legend = c('observed', 'polynomial fit'), lty = c(1, 2), bty = 'n')
    par(mfrow = c(1, 1))
  }
  return(list(fit = fit,
              sil.underwood = sil.underwood,
              sil.HomoTherm = sil.endoR,
              AREA = AREA))
}
