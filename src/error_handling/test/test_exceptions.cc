#include <UnitTest++.h>

#include "../exceptions.hh"

using namespace Exceptions;

void raise_Jali_exception ()
{

    Jali_exception exception;
    Jali_throw (exception);

}

SUITE (Exceptions)
{

    TEST (Behavior)
    {

        // Check the default value of RAISE
        CHECK_EQUAL (exception_behavior (), RAISE);

        set_exception_behavior_abort ();
        CHECK_EQUAL (exception_behavior (), ABORT);

        set_exception_behavior_raise ();
        CHECK_EQUAL (exception_behavior (), RAISE);

        set_exception_behavior (ABORT);
        CHECK_EQUAL (exception_behavior (), ABORT);

        set_exception_behavior (RAISE);
        CHECK_EQUAL (exception_behavior (), RAISE);

    }

    TEST (Catch)
    {

        CHECK_THROW (raise_Jali_exception (), Jali_exception);

    }


}
