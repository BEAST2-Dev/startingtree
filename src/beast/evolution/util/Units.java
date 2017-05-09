
package beast.evolution.util;


import beast.evolution.tree.TraitSet;

import java.io.Serializable;

/**
 * interface holding unit constants
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public interface Units extends Serializable {

    public enum Type {
//        SUBSTITUTIONS(XMLUnits.SUBSTITUTIONS), GENERATIONS(XMLUnits.GENERATIONS),
        DAY(TraitSet.Units.day.toString()), MONTH(TraitSet.Units.month.toString()),
        YEAR(TraitSet.Units.year.toString());

        Type(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }

        private final String name;
    }

    /**
     * @return the units for this object.
     */
    Type getUnits();

    /**
     * Sets the units for this object.
     *
     * @param units to use
     */
    void setUnits(Type units);


}